#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                                           !
#   ! Copyright 2021 by Tyler C. Sterling and Dmitry Reznik,                    !
#   ! University of Colorado Boulder                                            !
#   !                                                                           !
#   ! This file is part of the pynamic-structure-factor (PSF) software.         !
#   ! PSF is free software: you can redistribute it and/or modify it under      !
#   ! the terms of the GNU General Public License as published by the           !
#   ! Free software Foundation, either version 3 of the License, or             !
#   ! (at your option) any later version. PSF is distributed in the hope        !
#   ! that it will be useful, but WITHOUT ANY WARRANTY; without even the        !
#   ! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  !
#   ! See the GNU General Public License for more details.                      !
#   !                                                                           !
#   ! A copy of the GNU General Public License should be available              !
#   ! alongside this source in a file named gpl-3.0.txt. If not see             !
#   ! <http://www.gnu.org/licenses/>.                                           !
#   !                                                                           !
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import multiprocessing as mp
import numpy as np
from timeit import default_timer as timer
from scipy.fft import fft, fftfreq
import mod_io 
from mod_utils import print_stdout, PSF_exception
import mod_xlengths

class sqw:

    """
    store the scattered intensity arrays and calculator. see the doc string at the end of this file 
    for more info. calculates sqw, timeavg scattering, and bragg scattering if requested. see the 
    docstring at the end of this file for more info.

    """

    # ---------------------------------------------------------------------------------------
    
    def __init__(self,invars,Qpoints):

        """
        setup freq. grid from MD params and initialize variables to hold SQW
        i tried padding the FFT w 0's, but spectral leakage was waaayyy worse
        """

        # effective number of steps, i.e. num in each block that is read and computed
        self.block_steps = (invars.total_steps//invars.stride)//invars.num_blocks 
        self.pos = np.zeros((self.block_steps,invars.num_atoms,3)) # time-steps, atoms, xyz
        self.atom_types = np.zeros((self.block_steps,invars.num_atoms)).astype(int) # see mod_io
        self.box_lengths = [0,0,0] # read from traj file
        self.num_blocks = len(invars.blocks) # see mod_invars
        self.counter = 1 

        # tools to fill xlengths array for ins/xray scattering
        self.xlengths = np.zeros((self.block_steps,invars.num_atoms)) # this is an array
        self.xlengths_tools = mod_xlengths.scattering_lengths(invars.num_types) # this is a class

        # "rescale" the intensity to make plotting easier.
        self.common_rescale = 1e6

        # only create sqw array if requested. saves time to not do this if only bragg and/or timeavg
        if invars.compute_sqw:

            message = 'calculating the dynamical intensity'
            print_stdout(message,msg_type='NOTE')

            # setup frequency grid for fft
            self.num_freq = self.block_steps
            self.dt_eff = invars.dt*invars.stride*1e-3 # in units of ps, input is fs
            self.meV = fftfreq(self.num_freq,self.dt_eff)*4.13567
            self.df = self.meV[1]-self.meV[0]
            self.max_freq = self.meV.max()
            self.sqw_norm = self.num_freq # change this to change how time FT is normalized

            # print the energy resolution, max, etc
            message = (f'max freq: {self.max_freq:2.3f} meV\n'
                    f' number of freq.: {self.num_freq}\n'
                    f' resolution: {self.df:2.3e} meV\n')
            print_stdout(message,msg_type='frequency Grid')

            message = ('the time Fourier transform is normalized so that the MEAN\n'
                        ' of the dynamical intensity over all energies is equal to the\n'
                        ' total intensity averaged over time')
            print_stdout(message,msg_type='NORMALIZATION')

            # the sqw array
            self.sqw = np.zeros((self.num_freq,Qpoints.total_Qsteps))        

        # only create bragg array if requested
        if invars.compute_bragg:
            message = 'calculating the Bragg intensity'
            print_stdout(message,msg_type='NOTE')
            self.bragg = np.zeros(Qpoints.total_Qsteps)

        # only create timeavg array if requested
        if invars.compute_timeavg:
            message = 'calculating the time-avgeraged intensity'
            print_stdout(message,msg_type='NOTE')
            self.timeavg = np.zeros(Qpoints.total_Qsteps)

    # ----------------------------------------------------------------------------------------

    def calculate(self,invars,Qpoints,lattice,traj_file):

        """
        calculate SQW from MD data
        """

        # this call does all of the 'work'
        self._loop_over_blocks(invars,Qpoints,lattice,traj_file)

        # average over the blocks
        if invars.compute_sqw:
            self.sqw = self.sqw/self.num_blocks 

        if invars.compute_bragg:
            self.bragg = self.bragg/self.num_blocks 

        if invars.compute_timeavg:
            self.timeavg = self.timeavg/self.num_blocks 

        # free up memory (probably is already done by garbage collection)
        del self.pos, self.atom_types, self.xlengths

    # =======================================================================================
    # ------------------------------ private methods ----------------------------------------
    # =======================================================================================

    def _loop_over_blocks(self,invars,Qpoints,lattice,traj_file):

        """
        contains outer loop over blocks

        info about scattering lengths: there should be 1 length per TYPE, in order
        of types. e.g. for 4 types = 1,2,3,4 there should be for lengths atom 1 : length 1,
        atom 2 : lenght 2, etc... i am also assuming that dump_modify sort id was used so
        that the order of atoms is the  same for each step. this can be changed easily if
        not the case using the atom_types variable, but that will slow down the calc a little.
        the b_array variable has shape [num_steps, num_atoms] to vectorize calculating the
        neutron weighted density-density correlation fn
        """

        for block_index in invars.blocks: # loop over blocks to 'ensemble' average

            # used below
            self.block_index = block_index

            # print progress and start timer
            start_time = timer()
            message = '\n............................................'
            print_stdout(message)
            message = f' now on block {self.counter} out of {self.num_blocks}'
            print_stdout(message,msg_type='NOTE')

            # get the positions from file
            traj_file.parse_trajectory(invars,self) 

            # check that the number of b's defined in input file are consistent with traj file
            if np.unique(self.atom_types[0,:]).shape[0] != invars.num_types:
                message = 'number of types in input file doesnt match simulation'
                raise PSF_exception(message)

            # look up ins scattering lengths OR parameters to compute xray form factors.
            self.xlengths_tools.map_types_to_data(invars,self)

            # box lengths read from traj file
            a = self.box_lengths[0]/invars.supercell[0] 
            b = self.box_lengths[1]/invars.supercell[1]
            c = self.box_lengths[2]/invars.supercell[2]

            # print box lengths read from traj file to compare to input file
            message = f'cell lengths from hdf5 file: {a:2.3f} {b:2.3f} {c:2.3f} Angstrom'
            print_stdout(message,msg_type='NOTE')

            # recall, only ortho lattice vectors used (for now)
            if invars.recalculate_cell_lengths: # optionally recalculates from avg in MD file
                lattice.lattice_vectors = np.array([[a,0,0],[0,b,0],[0,0,c]])
                lattice.recompute_lattice()         # recompute reciprocal lattice
                Qpoints.reconvert_Q_points(lattice) # convert Q to 1/A in new basis

            Q_start_time = timer() # track time per Q not including the read/write time

            message = ('printing progess for process 0, which has >= the number of Q on other' \
                    ' processes.\n -- now entering loop over Q -- ')
            print_stdout(message,msg_type='NOTE')

            # -----------------------------------------------------------------------------------------
            #                 ------------- multiprocessing part. -------------
            # -----------------------------------------------------------------------------------------
            
            # a Queue to hold the retured SQW data
            self.mp_queue = mp.Queue()

            # a container to hold the processes
            procs = []

            # loop over processes, setting up the loop over Q on each.
            for pp in range(invars.num_processes):
                procs.append(mp.Process(target=self._loop_over_Q, args=(pp,invars,Qpoints)))

            # now start running the function on each process
            for proc in procs:
                proc.start()

            # note, doing it this way with the queue 'blocks' until the next processes adds to queue if 
            # it is empty. I dont know if this will freeze the whole calculation or just the background 
            # proc that is running the queue. anyway, everything with the queue has to be done before
            # joining the procs or the data will corrupt/crash the program.

            # get the stuff calculated on each proc
            for pp in range(invars.num_processes):

                sqw_pp, bragg_pp, timeavg_pp, proc = self.mp_queue.get()
                Q_inds = Qpoints.Q_on_procs[proc]

                # now put it into main arrays if requested
                if invars.compute_sqw:
                    self.sqw[:,Q_inds] = self.sqw[:,Q_inds]+sqw_pp
                if invars.compute_bragg:
                    self.bragg[Q_inds] = self.bragg[Q_inds]+bragg_pp
                if invars.compute_timeavg:
                    self.timeavg[Q_inds] = self.timeavg[Q_inds]+timeavg_pp

            # now close the queue and rejoin its proc
            self.mp_queue.close()
            self.mp_queue.join_thread()

            # wait here for all to finish before moving on to the next block
            for proc in procs:
                proc.join()

            # -----------------------------------------------------------------------------------------
            #                 ------------ end of multiprocessing part ----------------
            # -----------------------------------------------------------------------------------------

            # optionally save progress
            if invars.save_progress:

                if self.counter != self.num_blocks: # dont write if this is the last block

                    if invars.compute_sqw:
                        f_name = invars.outfile_prefix+f'_SQW_B{block_index}.hdf5'
                        mod_io.save_sqw(invars,Qpoints.reduced_Q,self.meV,self.sqw/self.counter,f_name)

                    if invars.compute_bragg:
                        f_name = invars.outfile_prefix+f'_BRAGG_B{block_index}.hdf5'
                        mod_io.save_bragg(invars,Qpoints.reduced_Q,self.bragg/self.counter,f_name)

                    if invars.compute_timeavg:
                        f_name = invars.outfile_prefix+f'_TIMEAVG_B{block_index}.hdf5'
                        mod_io.save_timeavg(invars,Qpoints.reduced_Q,self.timeavg/self.counter,f_name)

            # print timing to screen
            end_time = timer()
            elapsed_time = end_time-start_time
            Q_time = end_time-Q_start_time
            io_time = elapsed_time-Q_time

            # time per Qpoint
            Q_time = Q_time/len(Qpoints.Q_on_procs[0]) # avg over all Q
            message = f' avg time per Q-point:      {Q_time:2.3f} seconds'
            print_stdout(message,msg_type='TIMING')

            # time spent in i/o
            message = f' total io time:             {io_time:2.3f} seconds'
            print_stdout(message)

            # total time for in this method
            message = (f' total time for this block: {elapsed_time:2.3f} seconds'
                        f' ({elapsed_time/60:2.3f} minutes)')
            print_stdout(message)

            # update the block counter
            self.counter = self.counter+1

    # -----------------------------------------------------------------------------------------------

    def _loop_over_Q(self,proc,invars,Qpoints):

        """
        this target is passed to each mp Process. it loops over the Q points assigned to each 
        process and puts the results into the Queue. a background process collects all of this and 
        assembles to total S(Q,w), S(Q), etc. arrays. could be easily modified to work in serial.
        """

        # get the inds of which Q to do
        Q_inds = Qpoints.Q_on_procs[proc]

        # get how many there are
        num_Q = len(Q_inds) 

        # return dummy (None) sqw array if not requested
        if invars.compute_sqw:
            sqw_pp = np.zeros((self.num_freq,num_Q))
        else:
            sqw_pp = None

        # return dummy (None) bragg array if not requested
        if invars.compute_bragg:
            bragg_pp = np.zeros(num_Q)
        else:
            bragg_pp = None

        # return dummy (None) timeavg array if not requested
        if invars.compute_timeavg:
            timeavg_pp = np.zeros(num_Q)
        else:
            timeavg_pp = None

        # loop over the Q points
        for qq in range(num_Q):

            # print status if on proc 0
            if proc == 0:
                message = f' now on Q-point {qq+1} out of {num_Q}'
                print_stdout(message)

            # the Qpoint to do
            Q_ind = Q_inds[qq]
            Q = Qpoints.total_Qpoints[Q_ind,:].reshape((1,3)) # 1/Angstrom
            self.Q_norm = np.sqrt(Q[0,0]**2+Q[0,1]**2+Q[0,2]**2) # |Q|

            # if xray, need to compute f(|Q|), which are placed in self.xlengths
            if invars.exp_type == 'xray':
                self.xlengths_tools.compute_xray_form_fact(self,invars)

            # space FT by vectorized Q.r dot products and sum over atoms. (tile prepends new axes)
            exp_iQr = np.tile(Q,reps=[self.block_steps,invars.num_atoms,1])*self.pos # Q.r
            exp_iQr = np.exp(1j*exp_iQr.sum(axis=2))*self.xlengths # sum over x, y, z
            exp_iQr = exp_iQr.sum(axis=1) # sum over atoms

            # compute bragg intensity = |<rho(Q,t)>|**2
            if invars.compute_bragg:
                bragg_pp[qq] = np.abs((exp_iQr).mean())**2/self.common_rescale

            # compute timeavg intensity = <|rho(Q,t)|**2>
            if invars.compute_timeavg:
                timeavg_pp[qq] = (np.abs(exp_iQr)**2).mean()/self.common_rescale

            # compute dynamical intensity = |rho(Q,w)|**2
            if invars.compute_sqw:
                sqw_pp[:,qq] = np.abs(fft(exp_iQr))**2/self.sqw_norm/self.common_rescale

        # put this stuff into the 'Queue' so that some other process can put it into the main arrays
        self.mp_queue.put([sqw_pp,bragg_pp,timeavg_pp,proc])

    # ------------------------------------------------------------------------------------------------

    """

    description:

    1) the outer loop over 'blocks' (_loop_over_block) computes S(Q,w) for chunks 
    of the data in the specified file. the calculation in each block is independent of
    the others. they are averaged at the end. we read in the positions in the block at 
    the start of the outer loop. we also set up the scattering lengths array. if data 
    were 'wrapped' during the MD simulation, use the flag unwrap_positions = 1 to 'unwrap' 
    them (i.e. unimpose minimum image convention). if we don't unwrap, atoms near the 
    boundary of the simulation box can jump discontinously to the other side and the 
    Debye-Waller factor is screwed up.

    2) in each block, split the set of Q points between the requested number of processes
    and loop over the subset of Q points on each. the calculation at a given Q is independent
    of the others. (note, we could parallelize over blocks too, but since the calc. is already 
    memory intensive, it 

    3) at each Q point, compute the neutron-weighted density, rho (i.e. space FT weighted by 
    scattering lengths). we need the cross-correlation function of rho and its complex conjugate
    (sortof). <rho(0),rho(t)> ==  F(Q,t). this is sometimes called the intermediate scattering 
    function. the dynamic scattering function is time-FT[F(Q,t)] = S(Q,w). since x(t) is (approximated
    as) classical, everything commutes. then we just use (sortof) the Wiener-Khinchin theorem 
    to go directly from rho(t) to S(Q,w) = |time-FT[rho](Q,-w)|**2. The -w comes from the positive 
    time in rho(-Q,t). see my notes in the doc directory. my code basically returns S(Q,-w), but the 
    +/- w components are nearly identical. for refs, see Dove: "Lattice Dynamics," 
    Allen: "Computer Simulation of Liquids," and Squires: "Theory of Thermal Neutron Scattering"

    addendum) it also now can compute time averaged intensity and bragg (also timeaveraged) to 
    analyze diffuse (total-bragg) intensity

    validation:

    i didnt have another MD post-processing code to compare the output to (if i did, i wouldn't 
    have written this...) but i computed S(Q,w) in tersoff-silicon using this code. i also computed 
    the phonons using phonopy and a hack-job of an interface to lammps. i then used the eigenvectors 
    from phonopy in SNAXS, which computes S(Q,w) from the harmonic phonon expansion and my results
    match SNAXS **VERY** well. so i assert that this method is valid and accurate.

    notes:

    - the SQW array is normalized so that the AVERAGE over all frequencies/energies is equal 
    to the time averaged diffuse scattering intensity. this is because 1) the bragg/timeavg 
    intensities are averaged as opposed to integrated and 2) to make it easier to handle the 
    dimensions: integrated over ENERGY doesn't have the same dimensions as integrating over
    FREQUENCY
    - the space FT is all-ready highly vectorized. could probably be improved
    using some fancier LAPACK or BLAS functions to do vector products, but i dont think we can FFT
    the Q transform since its not on reduced q grid. maybe?
    
    this seems to run much faster than how i was doing with mpi4py. i think that was slow
    because each proc. was reading from the file concurrently. this way, we just read once
    and broadcast. I think this could be made to work with mpi4py too, but my hope with
    multiprocessing is to put the pos array into shared memory to reduce the memory
    footprint if possible
            
    - note on multiprocessing: doing it this way with the queue 'blocks' until the next processes 
    adds to queue if it is empty. I dont know if this will freeze the whole calculation or just 
    the background proc that is running the queue. anyway, everything with the queue has to be 
    done before joining the procs or the data will corrupt/crash the program.

    """

    # ---------------------------------------------------------------------------------------------




