import numpy as np
from timeit import default_timer as timer

class calc:

    def __init__(self,params):

        print(f'\n\tComputing S(Q,w) now. Check {params.log_file} for progress.\n')

        self.sqw = np.zeros((params.num_freq,params.Qsteps))
        self.pos = np.zeros((params.block_steps,params.num_atoms,3)) # time, atoms, xyz
        self.atom_ids = np.zeros((params.block_steps,params.num_atoms))
        self.b_array = np.zeros((params.block_steps,params.num_atoms))
        self.box_lengths = [-1,0,0]

        # this computes SQE
        self._loop_over_blocks(params)
        
        # should close hdf5's when done with them ... 
        params.traj_handle.close()

    def _loop_over_blocks(self,params):
        
        if params.debug == True:
            blocks = 1
        else:
            blocks = params.blocks
        for b in range(blocks):

            start_time = timer() # start the timer

            self.block_index = b
            
            # this part gets the positions and assembles the b_arrays
            # the optionally recomputes the lattice (and reciprocal lattice) vectors 
            # before computing the Fourier transforms
            if params.file_format == 'hdf5': # only support for hdf5 remains
                self._parse_traj_file(params)
                self._make_b_array(params)
                a = self.box_lengths[0]/params.supercell[0]
                b = self.box_lengths[1]/params.supercell[1]
                c = self.box_lengths[2]/params.supercell[2]

                params.log_handle.write(f'\n** NOTE: Box lengths from hdf5 file: {a:2.2f} {b:2.2f} {c:2.2f} Angstr **\n')
                params.log_handle.write('\n** NOTE: Box lengths from Input: {:2.2f} {:2.2f} {:2.2f} Angstr **\n'.format(
                    params.lattice_vectors[0][0],params.lattice_vectors[1][1],params.lattice_vectors[2][2]))

                if params.recalculate_box_lengths == True: # recalculates from avg in MD file
                    params.log_handle.write('\n** NOTE: Using box lengths from hdf5 file **\n')
                    params.log_handle.flush()
                    params.lattice_vectors = [[a,0,0],[0,b,0],[0,0,c]]
                    if params.Qpoints_file == False:
                        params._make_Qparams()
                    else:
                        params._Qpoints_from_list()
                else:
                    params.log_handle.write('\n** NOTE: Using box lengths from input **\n')
                    params.log_handle.flush()

            else:
                params.log_handle.write('\n** ERROR: Unkown file format "{}" **\n'.format(
                    params.file_format))
                params.log_handle.flush()
                exit()

            ### the actual calculation starts here
            self._loop_over_Q(params)

            end_time = timer()  # stop the timer
            elapsed_time = (end_time-start_time)/60 # minutes
            params.log_handle.write(f'\n Elapsed time for this block: '
                                    f'{elapsed_time:2.3f} minutes')
            params.log_handle.write(f'\n Time per Q-point: '
                                    f'{elapsed_time*60/params.Qsteps:2.3f} seconds\n')
            params.log_handle.flush()

        self.sqw = self.sqw/params.blocks # average it

    def _loop_over_Q(self,params):

        params.log_handle.write('\n** Loop over Q **\n')
        params.log_handle.flush()

        for q in range(params.Qsteps):
            self.q = q
            params.log_handle.write(f' Now on Qpoint {self.q+1} out of {params.Qsteps}\n')
            params.log_handle.flush()
            self.Q = params.Qpoints[q,:]
            self._compute_rho(params)

    def _compute_rho(self,params):

        exp_Qr = np.tile(self.Q.reshape(1,3),
            reps=[params.block_steps,params.num_atoms,1])*self.pos
        exp_Qr = np.exp(1j*exp_Qr.sum(axis=2))*self.b_array 
        self.sqw[:,self.q] = self.sqw[:,self.q]+np.abs(np.fft.fft(exp_Qr.sum(axis=1),
          n=params.num_freq))**2

        ### i think this is right. if we take the positions to be classical, then they commute and 
        ### everything can be factored willy nilly. then the S(Q,w) really can be written as a 
        ### product of the Fourier transforms of the neutron-weighted densities. 


    ################################################################################
    ################### auxiallary stuff below here ################################
    ################################################################################

    ## methods to parse pos files
    # i deleted the methods to read from text files. only use hdf5.

    def _parse_traj_file(self,params):

        inds = [self.block_index*params.block_steps,
                        (self.block_index+1)*params.block_steps]

        params.log_handle.write(f'\n** Reading Velocities ({params.file_format}) **\n')
        params.log_handle.write(f' Now on block {self.block_index+1} '
                        f'out of {params.blocks}\n')
        params.log_handle.write('\n Now reading:\n')
        params.log_handle.flush()

        self.box_lengths = [0,0,0]
        self.box_lengths[0] = np.mean(params.traj_handle['box_bounds'][inds[0]:inds[1],1]-
                params.traj_handle['box_bounds'][inds[0]:inds[1],0])
        self.box_lengths[1] = np.mean(params.traj_handle['box_bounds'][inds[0]:inds[1],3]-
                params.traj_handle['box_bounds'][inds[0]:inds[1],2])
        self.box_lengths[2] = np.mean(params.traj_handle['box_bounds'][inds[0]:inds[1],5]-
                params.traj_handle['box_bounds'][inds[0]:inds[1],4])
        params.log_handle.write(f'\tBox lengths: {self.box_lengths[0]:2.3f} '
            f'{self.box_lengths[1]:2.3f} {self.box_lengths[2]:2.3f} Angst.\n')
        params.log_handle.flush()

        self.pos[:,:,:] = params.traj_handle['pos_data'][inds[0]:inds[1],:,:]
        self.atom_ids[0,:] = params.traj_handle['atom_types'][:]

        if params.unwrap_pos == True:
            params.log_handle.write('\n** NOTE: Unwrapping positions **\n ')
            params.log_handle.flush()
            self._unwrap_positions(params)


    #################################################################
    # some other little methods to help with stuff above

    def _make_b_array(self,params):

        for a in range(params.num_atoms):
            self.b_array[0,a] = params.b[f'{self.atom_ids[0,a]:g}']
        self.b_array = np.tile(self.b_array[0,:].reshape(1,params.num_atoms),
                reps=[params.block_steps,1])


    ##################################################################

    def _unwrap_positions(self,params):
        
        # un-apply minimum image convention so that positions dont jump
        # discontinuously by a full box length
        # this wont be memory effecient but should be fast... ish

        lx = self.box_lengths[0]
        ly = self.box_lengths[1]
        lz = self.box_lengths[2]

        self.shift = self.pos-np.tile(self.pos[0,:,:].reshape(
            (1,params.num_atoms,3)),reps=[params.block_steps,1,1])

        dr = -lx*(self.shift[:,:,0] > lx/2).astype(int)
        dr = dr+lx*(self.shift[:,:,0] <= -lx/2).astype(int)
        self.shift[:,:,0] = np.copy(dr)

        dr = -ly*(self.shift[:,:,1] > ly/2).astype(int)
        dr = dr+ly*(self.shift[:,:,1] <= -ly/2).astype(int)
        self.shift[:,:,1] = np.copy(dr)

        dr = -lz*(self.shift[:,:,2] > lz/2).astype(int)
        dr = dr+lz*(self.shift[:,:,2] <= -lz/2).astype(int)
        self.shift[:,:,2] = np.copy(dr)
        
        self.pos = self.pos+self.shift
        
        del self.shift






