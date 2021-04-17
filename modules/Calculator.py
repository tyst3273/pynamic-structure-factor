
import numpy as np
from timeit import default_timer as timer
from FileIO import save_sqw, print_stdout

class calc:


    def __init__(self,params):

        """
        this class contains the data structures and methods relevant to SQW. 
        """

        self.sqw = np.zeros((params.num_freq,params.Qsteps))            # SQW array
        self.pos = np.zeros((params.block_steps,params.num_atoms,3))    # time-steps, atoms, xyz
        self.atom_ids = np.zeros((params.block_steps,params.num_atoms)) # used to map to scatt. lenghts
        self.b_array = np.zeros((params.block_steps,params.num_atoms))  # scattering lengths
        self.box_lengths = [0,0,0]                                      # read from traj file
        self.num_blocks = len(params.blocks)
        self.counter = 1


    # =============================================================
    # ****************    public methods    ***********************
    # =============================================================

    def run(self,params):

        """
        description:

        this a wrapper calling private methods to calculate SQW. 
        i put into a method and not in __init__() to make it eaier
        to write different methods later and just call those instead 


        an outline is as follows:

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

        # ================     outer loop over blocks     ===================

        self._loop_over_blocks(params)

        self.sqw = self.sqw/self.num_blocks # average over the blocks

        f_name = params.outfile_prefix+f'_P{params.my_rank}_BX.dat' # final file
        save_sqw(params,self.sqw,f_name=f_name)

        self._clean_up()
        params.clean_up()



    # ====================================================================
    # ********   primary private methods used to Compute SQW    **********
    # ====================================================================
    
    def _loop_over_blocks(self,params):

        """
        contains outer loop over blocks
        """

        for block_index in params.blocks:      

            if params.my_rank == 0:
                start_time = timer()        
                message = f'now on block {self.counter} out of {self.num_blocks}'
                print_stdout(message,msg_type='NOTE')

            self.block_index = block_index      
            
            # ================================================================
            # this part gets the positions and assembles the b_arrays
            # optionally recomputes the lattice and reciprocal-lattice vectors
            # and Q-points using the new reciprocal lattice vectors before 
            # computing the Fourier transforms. 
            # ================================================================

            self._parse_traj_file(params)
            self._make_b_array(params)
            a = self.box_lengths[0]/params.supercell[0] # read in the box lengths 
            b = self.box_lengths[1]/params.supercell[1]
            c = self.box_lengths[2]/params.supercell[2]

            # ==========   compute cell lengths and compare to input   ===============
            
            if params.my_rank == 0:
                message = f'cell lengths from hdf5 file: {a:2.3f} {b:2.3f} {c:2.3f} Angstrom'
                print_stdout(message,msg_type='NOTE')

            if params.recalculate_cell_lengths == True: # optionally recalculates from avg in MD file
                params.lattice_vectors = [[a,0,0],[0,b,0],[0,0,c]]
                params.convert_Q_to_1_over_A()          


            # ===============   loop over the Q points    ====================

            self._loop_over_Q(params) # sqw is calculated inside here. 


            # ==============    print timing to the log file     =============

            if params.my_rank == 0:
                end_time = timer()  
                elapsed_time = (end_time-start_time)/60 # minutes
                message = f' elapsed time for this block: {elapsed_time:2.3f} minutes'
                print_stdout(message,msg_type='TIMING')
                message = f' time per Q-point: {elapsed_time*60/params.Qsteps:2.3f} seconds'
                print_stdout(message)

            # =================================================================
            # sqw is added to each block. save each block, but divide by count  
            # when saving to average the current data
            # =================================================================

            if self.counter != self.num_blocks:

                f_name = params.outfile_prefix+f'_P{params.my_rank}_B{block_index}.dat'
                save_sqw(params,self.sqw/self.counter,f_name=f_name)

            self.counter = self.counter+1 # update the counter




    def _loop_over_Q(self,params):

        """
        contains inner loop over Q points
        """

        if params.my_rank == 0:
            message = ('printing progess for rank 0, which has >= the number of Q on other procs.\n'
                        ' -- now entering loop over Q -- ')
            print_stdout(message,msg_type='NOTE')

        for q in range(params.Qsteps):  # do the loop over Q points. 
            self.q = q
            self.Q = params.Qpoints[q,:]

            if params.my_rank == 0:
                message = f' now on Q-point {q+1} out of {params.Qsteps}'
                print_stdout(message)

            # ======   compute the neutron weighted density and correlation functions   ========

            self._compute_rho(params)




    def _compute_rho(self,params):

        """
        space-FT to get rho and use convolution theorem to compute 
        S(Q,w) = time-FT[<rho(t),rho(0)>] = |time-FT[rho(t)]|**2
        """

        # ======================================================================================
        # i think this is right. if we take the positions to be classical, then they commute and
        # everything can be factored willy nilly. then the S(Q,w) can be written as a product of
        # the time-FTs of the neutron-weighted densities (which are the space FTs of position).
        # ======================================================================================

        exp_iQr = np.tile(self.Q.reshape(1,3),
            reps=[params.block_steps,params.num_atoms,1])*self.pos
        exp_iQr = np.exp(1j*exp_iQr.sum(axis=2))*self.b_array 
        self.sqw[:,self.q] = self.sqw[:,self.q]+np.abs(np.fft.fft(exp_iQr.sum(axis=1),
          n=params.num_freq))**2




    # ====================================================================
    # ************  auxilliary private methods used above   **************
    # ====================================================================

    def _clean_up(self):

        del self.pos, self.b_array, self.atom_ids, 



    def _parse_traj_file(self,params):

        """
        this probably should have been put in the FileIO module

        parse the hdf5 file for positions, box lenghts, and atom types. 
        assemble b_array (scattering lengths) and  optionally call method
        to unwrap positions.
        """

        inds = [self.block_index*params.block_steps,    # list of the indicies in the block
                        (self.block_index+1)*params.block_steps]

        if params.my_rank == 0:
            message = 'now reading positions:'
            print_stdout(message,msg_type='NOTE')

        self.box_lengths = [0,0,0]      # get the box lengths
        self.box_lengths[0] = np.mean(params.traj_handle['box_bounds'][inds[0]:inds[1],1]-
                params.traj_handle['box_bounds'][inds[0]:inds[1],0])
        self.box_lengths[1] = np.mean(params.traj_handle['box_bounds'][inds[0]:inds[1],3]-
                params.traj_handle['box_bounds'][inds[0]:inds[1],2])
        self.box_lengths[2] = np.mean(params.traj_handle['box_bounds'][inds[0]:inds[1],5]-
                params.traj_handle['box_bounds'][inds[0]:inds[1],4])

        self.pos[:,:,:] = params.traj_handle['pos_data'][inds[0]:inds[1],:,:]   # get the positins 
        self.atom_ids[0,:] = params.traj_handle['atom_types'][:]                # get the atom TYPES


        # ======   (optionally) unimpose minimum image convention  =======

        if params.unwrap_pos == True:   

            if params.my_rank == 0:
                message = 'unwrapping positions'
                print_stdout(message,msg_type='NOTE')
            
            self._unwrap_positions(params)




    def _make_b_array(self,params):

        """
        assemble b_array: create a num_steps_in_block*num_atoms array with the scattering
        lengths mapped from b in params to atom_types read from hdf5
        """

        for a in range(params.num_atoms):
            self.b_array[0,a] = params.b[f'{self.atom_ids[0,a]:g}']
        self.b_array = np.tile(self.b_array[0,:].reshape(1,params.num_atoms),
                reps=[params.block_steps,1])




    def _unwrap_positions(self,params):

        """    
        un-apply minimum image convention so that positions dont jump discontinuously 
        by a full box length. this wont be memory effecient but should be fast... ish

        the issue is that, if the positions are written out 'wrapped' (e.g. lammps 
        x y z instead of xu yu zu), then atoms that cross the box boundary are wrapped
        back to the other side of the box. in that case, atoms near the box boundary will 
        occasionally jump discontinously by a full-box lenght. i did some checking and this
        screws up the debye waller factor at low Q (i.e. wave-length ~ the box). it was 
        only a minor effect, but this takes care of it and isn't super expensive. 

        the solution is to 'un-impose' the minimum image convention. i.e. treat every 
        atom as the center of the cell at t=0 and at all other times t', if the atom has 
        deviated by half a box length (i.e. outside the cell since the atom is at the 
        center), shift it back. See e.g. Allen: "Computer Simulation of Liquids".
        """

        lx = self.box_lengths[0] # average box lengths in block.
        ly = self.box_lengths[1]
        lz = self.box_lengths[2]

        # ======   build an array to be added to the positions to shift them  =======
        shift = self.pos-np.tile(self.pos[0,:,:].reshape((1,params.num_atoms,3)),
                reps=[params.block_steps,1,1])

        # ==========   check whether to shift atoms in the x direction     ==========
        dr = -lx*(shift[:,:,0] > lx/2).astype(int)
        dr = dr+lx*(shift[:,:,0] <= -lx/2).astype(int)
        shift[:,:,0] = np.copy(dr)

        # ==========   check whether to shift atoms in the y direction     ==========
        dr = -ly*(shift[:,:,1] > ly/2).astype(int)
        dr = dr+ly*(shift[:,:,1] <= -ly/2).astype(int)
        shift[:,:,1] = np.copy(dr)

        # ==========   check whether to shift atoms in the z direction     ==========
        dr = -lz*(shift[:,:,2] > lz/2).astype(int)
        dr = dr+lz*(shift[:,:,2] <= -lz/2).astype(int)
        shift[:,:,2] = np.copy(dr)
        
        # =====================     apply the shift    ==============================
        self.pos = self.pos+shift
        

















