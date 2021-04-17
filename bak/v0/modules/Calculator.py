
import numpy as np
from timeit import default_timer as timer

class calc:

    def __init__(self,params):

        """
        this class contains the data structures and methods relevant to SQW. 
        """

        print(f'\n\tComputing S(Q,w) now. Check {params.log_file} for progress.\n')

        self.sqw = np.zeros((params.num_freq,params.Qsteps))            # SQW array
        self.pos = np.zeros((params.block_steps,params.num_atoms,3))    # time-steps, atoms, xyz
        self.atom_ids = np.zeros((params.block_steps,params.num_atoms)) # used to map to scattering lenghts
        self.b_array = np.zeros((params.block_steps,params.num_atoms))  # scattering lengths
        self.box_lengths = [-1,0,0]                                     # used if calculating cell lengths from file
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
        tersoff-silicon using this code. i also computed the phonons using phonopy and 
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

        self.sqw = self.sqw/self.num_blocks   # average over the blocks




    # ====================================================================
    # ********   primary private methods used to Compute SQW    **********
    # ====================================================================
    
    def _loop_over_blocks(self,params):

        """
        contains outer loop over blocks
        """

        for block_index in params.blocks:      
            start_time = timer()                # start a timer to print wall time to log file
            self.block_index = block_index      
            
            # ================================================================
            # this part gets the positions and assembles the b_arrays
            # optionally recomputes the lattice and reciprocal-lattice vectors
            # and Q-points using the new reciprocal lattice vectors before 
            # computing the Fourier transforms. 
            # ================================================================

            if params.file_format == 'hdf5':    # only support for hdf5 remains
                self._parse_traj_file(params)
                self._make_b_array(params)
                a = self.box_lengths[0]/params.supercell[0] # read in the box lengths 
                b = self.box_lengths[1]/params.supercell[1]
                c = self.box_lengths[2]/params.supercell[2]

                # ==========     compute cell lengths and compare to input    ===============

                params.log_handle.write(f'\n** NOTE: cell lengths from hdf5 file: {a:2.2f} {b:2.2f} {c:2.2f} Angstr **\n')
                params.log_handle.write('\n** NOTE: cell lengths from input: {:2.2f} {:2.2f} {:2.2f} Angstr **\n'.format(
                    params.lattice_vectors[0][0],params.lattice_vectors[1][1],params.lattice_vectors[2][2]))

                if params.recalculate_cell_lengths == True: # optionally recalculates from avg in MD file
                    params.log_handle.write('\n** NOTE: Using cell lengths from hdf5 file **\n')
                    params.lattice_vectors = [[a,0,0],[0,b,0],[0,0,c]]
                    params.gen_Qpoints()                    # compte Q points with these cell lengths
                else:                                       # otherwise just print to log
                    params.log_handle.write('\n** NOTE: Using cell lengths from input **\n')
                params.log_handle.flush()

            else:   # this useless since hdf5 is hard-enforced in the params.__init__() but left here for legacy reasons
                params.log_handle.write('\n** ERROR: Unkown file format "{}" **\n'.format(params.file_format))
                params.log_handle.flush()
                exit()


            # ===============   loop over the Q points    ====================

            self._loop_over_Q(params) # sqw is calculated inside here. 


            # ==============    print timing to the log file     =============

            end_time = timer()  # stop the timer
            elapsed_time = (end_time-start_time)/60 # minutes
            params.log_handle.write(f'\n Elapsed time for this block: '
                                    f'{elapsed_time:2.3f} minutes')
            params.log_handle.write(f'\n Time per Q-point: '
                                    f'{elapsed_time*60/params.Qsteps:2.3f} seconds\n')
            params.log_handle.flush()

            self.counter = self.counter+1



    def _loop_over_Q(self,params):

        """
        contains inner loop over Q points
        """

        params.log_handle.write('\n** Loop over Q **\n')
        params.log_handle.flush()

        for q in range(params.Qsteps):  # do the loop over Q points. 
            self.q = q
            params.log_handle.write(f' Now on Qpoint {self.q+1} out of {params.Qsteps}\n')
            params.log_handle.flush()   # print progress to log file. useful to check on long runs. 
            self.Q = params.Qpoints[q,:]


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

    def _parse_traj_file(self,params):

        """
        parse the hdf5 file for positions, box lenghts, and atom types. 
        assemble b_array (scattering lengths) and  optionally call method
        to unwrap positions.
        """

        inds = [self.block_index*params.block_steps,    # list of the indicies in the block
                        (self.block_index+1)*params.block_steps]

        params.log_handle.write(f'\n** Reading Velocities ({params.file_format}) **\n')
        params.log_handle.write(f' Now on block {self.counter} out of {self.num_blocks}\n')
        params.log_handle.write('\n Now reading:\n')
        params.log_handle.flush()

        self.box_lengths = [0,0,0]      # get the box lengths
        self.box_lengths[0] = np.mean(params.traj_handle['box_bounds'][inds[0]:inds[1],1]-
                params.traj_handle['box_bounds'][inds[0]:inds[1],0])
        self.box_lengths[1] = np.mean(params.traj_handle['box_bounds'][inds[0]:inds[1],3]-
                params.traj_handle['box_bounds'][inds[0]:inds[1],2])
        self.box_lengths[2] = np.mean(params.traj_handle['box_bounds'][inds[0]:inds[1],5]-
                params.traj_handle['box_bounds'][inds[0]:inds[1],4])
        params.log_handle.write(f' Box lengths: {self.box_lengths[0]:2.3f} '
            f'{self.box_lengths[1]:2.3f} {self.box_lengths[2]:2.3f} Angst.\n')
        params.log_handle.flush()

        self.pos[:,:,:] = params.traj_handle['pos_data'][inds[0]:inds[1],:,:]   # get the positins 
        self.atom_ids[0,:] = params.traj_handle['atom_types'][:]                # get the atom TYPES


        # ======   (optionally) unimpose minimum image convention  =======

        if params.unwrap_pos == True:   
            params.log_handle.write('\n** NOTE: Unwrapping positions **\n ')
            params.log_handle.flush()
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
        self.shift = self.pos-np.tile(self.pos[0,:,:].reshape((1,params.num_atoms,3)),
                reps=[params.block_steps,1,1])

        # ==========   check whether to shift atoms in the x direction     ==========
        dr = -lx*(self.shift[:,:,0] > lx/2).astype(int)
        dr = dr+lx*(self.shift[:,:,0] <= -lx/2).astype(int)
        self.shift[:,:,0] = np.copy(dr)

        # ==========   check whether to shift atoms in the y direction     ==========
        dr = -ly*(self.shift[:,:,1] > ly/2).astype(int)
        dr = dr+ly*(self.shift[:,:,1] <= -ly/2).astype(int)
        self.shift[:,:,1] = np.copy(dr)

        # ==========   check whether to shift atoms in the z direction     ==========
        dr = -lz*(self.shift[:,:,2] > lz/2).astype(int)
        dr = dr+lz*(self.shift[:,:,2] <= -lz/2).astype(int)
        self.shift[:,:,2] = np.copy(dr)
        
        # =====================     apply the shift    ==============================
        self.pos = self.pos+self.shift
        
        del self.shift # its a fairly large array; free up the memory.   






