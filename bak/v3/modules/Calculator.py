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
            # the hdf5 case recomputes the lattice (and reciprocal lattice) vectors 
            # before computing the Fourier transforms
            if params.file_format == 'default':
                self._parse_traj_file(params)
                self._make_b_array(params)
            elif params.file_format == 'ertekin':
                self._parse_traj_file_ertekin(params)
                self._make_b_array_ertekin(params)
            elif params.file_format == 'hdf5':
                self._parse_traj_file_hdf5(params)
                self._make_b_array(params)
                a = self.box_lengths[0]/params.supercell[0]
                b = self.box_lengths[1]/params.supercell[1]
                c = self.box_lengths[2]/params.supercell[2]
                params.lattice_vectors = [[a,0,0],[0,b,0],[0,0,c]]
                params._make_Qparams()
            else:
                params.log_handle.write('\n** ERROR: Unkown file format "{}" **\n'.format(
                    params.file_format))
                params.log_handle.flush()
                exit()

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

        # since we are leaving the 'phase factor' exp(-iwt) out in E.21 (see the comment
        # in the _parse_traj_file method), the complex conjugte of rho really is 
        # rho(-Q,t), so we can just square the FT of rho(Q)
        
        exp_Qr = np.tile(self.Q.reshape(1,3),
            reps=[params.block_steps,params.num_atoms,1])*self.pos
        exp_Qr = np.exp(1j*exp_Qr.sum(axis=2))*self.b_array 
        self.sqw[:,self.q] = self.sqw[:,self.q]+np.abs(np.fft.fft(exp_Qr.sum(axis=1),
          n=params.num_freq)*np.fft.fft(np.conj(exp_Qr).sum(axis=1),n=params.num_freq))**2

        ### Note: this might not be right. Before, I used the convolution theorem which 
        # was dubious... I double checked and the transform Kernel isn't conjugated too. So 
        # what is above now is actually better, except the output is complex with negative
        # values in the real part. The norm squared looks right and is actually better than
        # the old way I was doing it in v1, probably just because I am squaring the
        # intensity. I will need to double check the theory.

        ### new note: just do the math.. anyway, the phase factor in E.21 times the kernel
        # is exp(-iw(t'-t'')). We can change variables and integrate over t''-t' and then I
        # think the convolution theorem is okay to use? 



    ### auxiallary stuff below here

    def _parse_traj_file(self,params):
        """
        Need to consider something: if the positions are 'wrapped' when written to the
        file, the the positions will hop discountinously across to the other side of the
        box. A way to fix this: take the positions at t=0 as the center of the coordinates 
        and apply minimum image on the positions at arbitrary t.

        As of now, my positions are NOT wrapped when written to the file.

        Also need to consider 'phase factors,' see Dove E.21. Without working out the math,
        my thought is this: in the integral in E.23 has rho*exp(-iwt'), where rho contains
        exp(-iwt). Then you get integral( sum_j bj exp(iQ.rj(t) exp(-iw(t+t') dt'}
        Letting t+t' = tau, dtau = dt' and we can just Fourier transform the sum minus the
        phase factor. Lets try that first.
        """

        self.box_lengths = [-1,0,0]
        params.log_handle.write('\n** Reading Velocities **\n')
        params.log_handle.write(f' Now on block {self.block_index+1} '
                        f'out of {params.blocks}\n')        
        params.log_handle.flush()

        if self.box_lengths[0] == -1:
            for i in range(5):
                params.traj_handle.readline()
            tmp = params.traj_handle.readline().strip().split()
            self.box_lengths[0] = float(tmp[1])-float(tmp[0])
            tmp = params.traj_handle.readline().strip().split()
            self.box_lengths[1] = float(tmp[1])-float(tmp[0])
            tmp = params.traj_handle.readline().strip().split()
            self.box_lengths[2] = float(tmp[1])-float(tmp[0])
            params.log_handle.write(f' Box lengths: {self.box_lengths[0]:2.3f} '
                        f'{self.box_lengths[1]:2.3f} {self.box_lengths[2]:2.3f} Angst.\n')
            params.log_handle.flush()
            params.traj_handle.readline()
            for i in range(params.num_atoms):
                tmp = params.traj_handle.readline().strip().split()
                self.atom_ids[0,i] = int(tmp[1])
            params.traj_handle.seek(0)

        params.log_handle.write('\n Now reading:\n')
        params.log_handle.flush()
        for step in range(params.block_steps):
            if step%1000 == 0:
                params.log_handle.write(f' Now on step {step} out '
                            f'of {params.block_steps}\n')
                params.log_handle.flush()
            for i in range(9): # header
                params.traj_handle.readline()
            for atom in range(params.num_atoms):
                self.pos[step,atom,:] = params.traj_handle.readline().strip().split()[2:]


    ##########################################

    def _parse_traj_file_ertekin(self,params):
        """
        i can just write an *unwrap_positions* method that will take the block of
        trajectory data and apply minimum image to each atom as a function of time. 
        if it move further than a box size in any step, just add the box to it. we can 
        ignore this for now
        """

        params.log_handle.write(f'\n** Reading Velocities ({params.file_format}) **\n')
        params.log_handle.write(f' Now on block {self.block_index+1} '
                        f'out of {params.blocks}\n')
        params.log_handle.flush()

        self.box_lengths = [0,0,0]
        params.log_handle.write('\n Now reading:\n')
        params.log_handle.flush()

        for step in range(params.block_steps):

            if step%500 == 0:
                params.log_handle.write(f' Now on step {step} out '
                            f'of {params.block_steps}\n')
                params.log_handle.flush()

                for i in range(5):
                    params.traj_handle.readline()
                tmp = params.traj_handle.readline().strip().split()
                self.box_lengths[0] = float(tmp[1])-float(tmp[0])
                tmp = params.traj_handle.readline().strip().split()
                self.box_lengths[1] = float(tmp[1])-float(tmp[0])
                tmp = params.traj_handle.readline().strip().split()
                self.box_lengths[2] = float(tmp[1])-float(tmp[0])
                params.traj_handle.readline()
                params.log_handle.write(f'\tBox lengths: {self.box_lengths[0]:2.3f} '
                    f'{self.box_lengths[1]:2.3f} {self.box_lengths[2]:2.3f} Angst.\n')
                params.log_handle.flush()

            else:
                for i in range(9): # header
                    params.traj_handle.readline()
            
            for atom in range(params.num_atoms):
                tmp = params.traj_handle.readline().strip().split()
                self.atom_ids[step,atom] = int(tmp[1])
                self.pos[step,atom,:] = tmp[2:5]


    ##########################################

    def _parse_traj_file_hdf5(self,params):
        """
        still haven't wrapped the positions ... maybe it isn't necessary
        """

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


    #######################################################
    ######################################################

    def _make_b_array(self,params):

        for a in range(params.num_atoms):
            self.b_array[0,a] = params.b[f'{self.atom_ids[0,a]:g}']
        self.b_array = np.tile(self.b_array[0,:].reshape(1,params.num_atoms),
                reps=[params.block_steps,1])

    def _make_b_array_ertekin(self,params):

        for i in range(params.block_steps):
            for a in range(params.num_atoms):
                self.b_array[i,a] = params.b[f'{self.atom_ids[i,a]:g}']






