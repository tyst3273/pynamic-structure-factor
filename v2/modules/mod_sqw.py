import numpy as np
from timeit import default_timer as timer
from scipy.fft import fft
import mod_io 
from mod_utils import print_stdout, PSF_exception

class sqw:

    """
    store the SQW data and calculator
    """

    # ---------------------------------------------------------------------------------------
    
    def __init__(self,invars,Qpoints,rank):
        """
        setup freq. grid from MD params and initialize variables to hold SQW
        """
        self.rank = rank

        # effective number of steps, i.e. num in each block that is read and computed
        self.block_steps = (invars.total_steps//invars.stride)//invars.num_blocks 
        self.max_freq = 1e-12/invars.dt/invars.stride*4.13567
        self.num_freq = self.block_steps
        self.meV = np.linspace(0,self.max_freq,self.num_freq) 
        self.df = self.meV[1]-self.meV[0] # freq. resolution

        if self.rank == 0:
            message = (f'max freq: {self.max_freq/2:2.3f} meV\n number of freq.: {self.num_freq}\n'
                       f' resolution: {self.df:2.3e} meV\n')
            print_stdout(message,msg_type='frequency Grid')

        self.sqw = np.zeros((self.num_freq,Qpoints.Qsteps))        # SQW array
        self.pos = np.zeros((self.block_steps,invars.num_atoms,3)) # time-steps, atoms, xyz
        self.atom_ids = np.zeros((self.block_steps,invars.num_atoms)).astype(int)
        self.b_array = np.zeros((self.block_steps,invars.num_atoms))  
        self.box_lengths = [0,0,0] # read from traj file
        self.num_blocks = len(invars.blocks)
        self.counter = 1

    # ----------------------------------------------------------------------------------------

    def calculate(self,invars,Qpoints,lattice,traj_file):
        """
        calculate SQW from MD data
        """
        self._loop_over_blocks(invars,Qpoints,lattice,traj_file)

        self.sqw = self.sqw/self.num_blocks # average over the blocks

        if invars.save_progress:
            f_name = invars.outfile_prefix+f'_P{self.rank}_BF.hdf5' # final file
            mod_io.save_sqw(invars,Qpoints.reduced_Q,self.meV,self.sqw,f_name)

        del self.pos, self.b_array, self.atom_ids

    # =======================================================================================
    # ------------------------------ private methods ----------------------------------------
    # =======================================================================================

    def _loop_over_blocks(self,invars,Qpoints,lattice,traj_file):
        """
        contains outer loop over blocks
        """
        for block_index in invars.blocks:

            self.block_index = block_index

            if self.rank == 0:
                start_time = timer()
                message = f'now on block {self.counter} out of {self.num_blocks}'
                print_stdout(message,msg_type='NOTE')

            # get the positions from file
            traj_file.parse_trajectory(invars,self) 

            # -------------------- set up array of scattering lengths -------------------------
            # note, these should be 1 length per TYPE, in order of types. e.g. for 4 types,
            # 1,2,3,4 there should be for lenghts  atom 1 : length 1, atom 2 : lenght 2, etc...
            # also assuming that dump_modify sort id was used so that the order of atoms is the 
            # same for each step. this can be changed easily if not the case. 
            # ---------------------------------------------------------------------------------
            for aa in range(invars.num_atoms):
                self.b_array[0,aa] = invars.b[self.atom_ids[0,aa]-1]
            self.b_array = np.tile(self.b_array[0,:].reshape(1,invars.num_atoms),
                    reps=[self.block_steps,1])

            # box lenghts read from traj file
            a = self.box_lengths[0]/invars.supercell[0] 
            b = self.box_lengths[1]/invars.supercell[1]
            c = self.box_lengths[2]/invars.supercell[2]

            if self.rank == 0:
                message = f'cell lengths from hdf5 file: {a:2.3f} {b:2.3f} {c:2.3f} Angstrom'
                print_stdout(message,msg_type='NOTE')

            if invars.recalculate_cell_lengths: # optionally recalculates from avg in MD file
                invars.lattice_vectors = np.array([[a,0,0],[0,b,0],[0,0,c]])
                lattice.recompute_lattice(invars)
                Qpoints.reconvert_Q_points(lattice)

            if self.rank == 0:
                message = ('printing progess for rank 0, which has >= the number of Q on other procs.\n'
                            ' -- now entering loop over Q -- ')
                print_stdout(message,msg_type='NOTE')

            # do the loop over Q points
            for qq in range(Qpoints.Qsteps):  

                if self.rank == 0:
                    message = f' now on Q-point {qq+1} out of {Qpoints.Qsteps}'
                    print_stdout(message)

                Q = Qpoints.Qpoints[qq,:] # these are the ones in 1/Angstrom

                exp_iQr = np.tile(Q.reshape(1,3),reps=[self.block_steps,invars.num_atoms,1])*self.pos
                exp_iQr = np.exp(1j*exp_iQr.sum(axis=2))*self.b_array
                self.sqw[:,qq] = self.sqw[:,qq]+np.abs(fft(exp_iQr.sum(axis=1)))**2

            # print timing to the log file
            if self.rank == 0:
                end_time = timer()

                elapsed_time = (end_time-start_time)/60 # minutes
                message = f' elapsed time for this block: {elapsed_time:2.3f} minutes'
                print_stdout(message,msg_type='TIMING')

                message = f' time per Q-point: {elapsed_time*60/Qpoints.Qsteps:2.3f} seconds'
                print_stdout(message)

            # optionally save progress
            if invars.save_progress:
                if self.counter != self.num_blocks:
                    f_name = invars.outfile_prefix+f'_P{self.rank}_B{block_index}.hdf5'
                    mod_io.save_sqw(invars,Qpoints.reduced_Q,self.meV,self.sqw/self.counter,f_name)

            self.counter = self.counter+1 # update the counter





        """
        description:

        1) the outer loop over 'blocks' computes SQW for chunks of the data
        in the specified file. the calculation in each block is independent of
        the others. they are averaged at the end. we read in the positions
        in the block at the start of the outer loop. we also set up the
        scattering lengths array.

        2) in each block, we loop over Q points. the calculation at a given
        Q is independent of the others.

        3) at each Q point, compute the neutron-weighted density, rho (i.e. space
        FT weighted by scattering lengths). we need its auto-correlation function:
        <rho(t),rho(o)> == S(Q,t) == F(Q,t). this is sometimes called the intermediate
        scattering function. the dynamic scattering function is time-FT[F(Q,t)] = S(Q,w).
        since x(t) is classical, everything commutes. then we just use the convolution
        theorem to go directly from rho(t) to S(Q,w) = |time-FT[rho]|**2

        for refs, see Dove: "Lattice Dynamics," Allen: "Computer Simulation of Liquids,"
        and Squires: "Theory of Thermal Neutron Scattering"

        validation:

        i didnt have another MD post-processing code to compare the output to
        (if i did, i wouldn't have written this...) but i computed S(Q,w) in
        tersoff-silicon using this code. i also computed the phonons using phonopy and a
        hack-job of an interface to lammps. i then used the eigenvectors from phonopy
        in SNAXS, which computes S(Q,w) from the harmonic phonon expansion and my results
        match SNAXS **VERY** well. so i assert that this method is valid and accurate.

        notes:

        1) could parallelize over the blocks.
        2) could parallelize over Q points.
        3) the space FT is all-ready highly vectorized. could probably be improved
        using some fancier LAPACK or BLAS functions to do vector products, but the bottle
        neck for now is serial execution over blocks and Q-points.

        """


