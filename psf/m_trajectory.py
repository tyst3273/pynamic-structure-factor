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

import numpy as np
from psf.m_error import crash


class c_trajectory:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm,timers):

        """
        class for handling MD trajectory I/O. reads different file types. it's modular
        but a little clumsy to add new parsers to.

        really, new parsers should be their own classes that return the pos and types
        arrays in the expected format. oh well. that can be implemented later!
        """

        msg = '\n*** NOTE ***\n'
        msg += 'consider using multiprocessing.shared_memory to place the trajectory into\n'
        msg += 'shared memory that can later be accessed by forked processes!\n'
        print(msg)

        # copy refs to stuff
        self.comm = comm
        self.config = config
        self.timers = timers

        # the trajectory file
        self.num_atoms = self.config.md_num_atoms
        self.trajectory_format = self.config.trajectory_format
        self.md_time_step = self.config.md_time_step    
        self.effective_time_step = self.config.trajectory_stride*self.md_time_step
        self.unwrap_trajectory = self.config.unwrap_trajectory

        # set up the 'blocks' of indices for calculating on
        self._get_block_inds()

        _eff_step = self.effective_time_step/1000 # ps
        _duration = self.num_block_steps*_eff_step
        msg = f'effective timestep: {_eff_step} (ps)\n'
        msg += f'trajectory duration: {_duration:6.2f} (ps)'
        print(msg)

        # initialize the positions; they are cartesian coords in whatever units are in the file
        self.pos = np.zeros((self.num_block_steps,self.num_atoms,3))

        # determine whether or not to read box vectors from file; not need if not unwrapping
        if not self.unwrap_trajectory:
            self.read_box = False
            self.box_vectors = None

        # if unwrapping, check wheter or not to get from file
        else:
            if self.config.box_vectors is None: # read from traj file
                print('\ngetting box vectors from file')
                self.read_box = True
                self.box_vectors = np.zeros((3,3),dtype=float)
            else: # use input arg
                self.read_box = False
                self.box_vectors = self.config.box_vectors

        # if setting types and pos externally, dont do anything else here
        if self.trajectory_format in ['external']:
            return

        # the file with traj in it
        self.trajectory_file = self.config.trajectory_file

        # get atom types once and for all (if possible)
        self.get_atom_types()

    # ----------------------------------------------------------------------------------------------

    def _get_block_inds(self):

        """
        setup indices for block averaging
        """

        # inds in file
        _inds = np.arange(self.config.trajectory_skip,
            self.config.md_num_steps-self.config.trajectory_trim,self.config.trajectory_stride)
        _num_steps = _inds.size

        # split up remaining inds
        _num_blocks = self.config.num_trajectory_blocks

        # averging is done in frequency space; having blocks of data that are on
        # different frequency grids will require interpolation to average together
        # i want to avoid doing that so i crop the dataset so that all blocks
        # have same frequency
        if _num_steps % _num_blocks != 0:
            msg = '\n*** warning! ***\n'
            msg += 'number of_steps must be divisible by num_trajectory_blocks\n'
            msg += 'discarding the extra timesteps\n'
            print(msg)

        self.blocks = self.config.trajectory_blocks

        # get the steps per block
        _rem = _num_steps % _num_blocks
        self.num_block_avg = self.config.num_block_avg
        self.num_steps_used = _num_steps-_rem
        self.num_block_steps = self.num_steps_used // _num_blocks

        # set up block inds
        self.block_inds = np.zeros((_num_blocks,self.num_block_steps),dtype=int)
        for bb in range(_num_blocks):
            self.block_inds[bb,:] = _inds[bb*self.num_block_steps:(bb+1)*self.num_block_steps]

        msg = f'\n*** blocks ***\nnum_block_avg: {self.num_block_avg}\n'
        msg += f'num_block_steps: {self.num_block_steps}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def get_atom_types(self):

        """
        get an array that holds types of atom. only do for 1st step ...

        it is MANDATORY that data are sorted in the exact same way each timestep... the
        code assumes the order of the types of atoms doesnt change over time, so only calculates
        the scattering meta-data (e.g. xray form factors) at the first timestep
        """

        if self.trajectory_format == 'lammps_hdf5':
            self._read_types_lammps_hdf5()
        if self.trajectory_format == 'user_hdf5':
            self._read_types_user_hdf5()

        # do generic error checking
        _types = np.unique(self.types)
        if _types.max() >= self.comm.config.num_types or _types.min() < 0 \
            or _types.size != self.comm.config.num_types:
            msg = 'types in trajectory file are incompatible with types in input file\n'
            crash(msg)

    # ----------------------------------------------------------------------------------------------

    def set_external_types(self,types):

        """
        get an array that holds types of atom. given explicitly, then checked
        """

        self.types = np.array(types)

        _types = np.unique(self.types)
        if _types.max() >= self.comm.config.num_types or _types.min() < 0 \
            or _types.size != self.comm.config.num_types:
            msg = 'given external types are incompatible with types in input file\n'
            crash(msg)

    # ----------------------------------------------------------------------------------------------

    def set_external_pos(self,pos):

        """
        get an array that holds time-series of positions. given explicitly, then checked
        """

        self.external_pos = np.array(pos)

        msg = 'external positions have wrong shape\n'
        _shape = self.external_pos.shape
        if len(_shape) == 2:
            self.external_pos.shape = [1,_shape[0],_shape[1]]

        if self.external_pos.shape[0] != self.config.md_num_steps:
            crash(msg)
        if self.external_pos.shape[1] != self.config.md_num_atoms:
            crash(msg)
        if self.external_pos.shape[2] != 3:
            crash(msg)

    # ----------------------------------------------------------------------------------------------

    def read_trajectory_block(self,block_index=0):

        """
        get a block of trajectory data from the file and preprocess it as much as possible.
        note: block_index is the only argument. it labels which 'block' of trajectory data to get.
        the indices of the steps in the block (accounting for stride, strip, and trim) are already
        calculated in _get_block_inds() above.
        """

        # step indices to get from file for this block
        inds = self.block_inds[block_index,:]

        # --- add new file type parsers here! ---
        if self.trajectory_format == 'lammps_hdf5':
            self._read_pos_lammps_hdf5(inds)            
        if self.trajectory_format == 'user_hdf5':
            self._read_pos_user_hdf5(inds)
        elif self.trajectory_format == 'external':
            self._read_pos_external(inds)
        # ---------------------------------------

        # check if simulation box is orthorhombic; only matters if unwrapping
        if self.unwrap_trajectory:

            self.box_is_ortho = True
            _tmp = np.abs(self.box_vectors)
            _tmp[0,0] = 0.0; _tmp[1,1] = 0.0; _tmp[2,2] = 0.0
            if np.any(_tmp > 1e-3):
                self.box_is_ortho = False
        
            if self.read_box:
                msg = '\n*** box vectors from file ***\n'
                for ii in range(3):
                    msg += '\n'
                    for jj in range(3):
                        msg += f'{self.box_vectors[ii,jj]: 10.6f} '            

            self._unwrap_positions()

    # ----------------------------------------------------------------------------------------------
    
    def _unwrap_positions(self):

        """
        un-apply minimum image convention so that positions dont jump discontinuously
        by a full box length. this wont be memory effecient but should be fast... ish

        see 'note about unwrapping' in the README
        """

        if not self.box_is_ortho:
            msg = 'unwrapping positions is not supported for non-orthorhombic simulation\n' \
                  'boxes yet. either use a different box or turn off unwrapping (okay..)\n'
            crash(msg)

        self.timers.start_timer('unwrap_positions',units='s')

        msg = '\nunwrapping positions!\n'
        print(msg)

        # get box lengths... assumes the simulation cell and underlying unitcell are orthorhombic
        lx = self.box_vectors[0,0]
        ly = self.box_vectors[1,1]
        lz = self.box_vectors[2,2]

        # build an array to be added to the positions to shift them 
        shift = self.pos-np.tile(self.pos[0,:,:].reshape((1,self.num_atoms,3)),
                reps=[self.num_block_steps,1,1])

        # check whether to shift atoms in the x direction     
        dr = -lx*(shift[:,:,0] > lx/2).astype(int)
        dr = dr+lx*(shift[:,:,0] <= -lx/2).astype(int)
        shift[:,:,0] = np.copy(dr)

        # check whether to shift atoms in the y direction  
        dr = -ly*(shift[:,:,1] > ly/2).astype(int)
        dr = dr+ly*(shift[:,:,1] <= -ly/2).astype(int)
        shift[:,:,1] = np.copy(dr)

        # check whether to shift atoms in the z direction   
        dr = -lz*(shift[:,:,2] > lz/2).astype(int)
        dr = dr+lz*(shift[:,:,2] <= -lz/2).astype(int)
        shift[:,:,2] = np.copy(dr)

        # apply the shift    
        self.pos = self.pos+shift

        self.timers.stop_timer('unwrap_positions')

    # ----------------------------------------------------------------------------------------------

    def _read_types_user_hdf5(self):

        """
        get from user created hdf5 file. 
        """

        self.types = np.zeros(self.num_atoms,dtype=int)

        try:
            import h5py
            with h5py.File(self.trajectory_file,'r') as in_db:

                # it is assumed that the types are consecutive integers starting at 1
                # so that 1 is subtracted to match the python index starting at 0 ...
                self.types[:] = in_db['types'][:]-1

            # error check
            if self.num_atoms != self.types.size:
                msg = 'number of atoms in file doesnt match the input info.'
                msg += 'correct the configuration or check your trajectory file!\n'
                crash(msg)

        except Exception as _ex:
            msg = f'the hdf5 file\n  \'{self.trajectory_file}\'\ncould not be read!\n'
            msg += 'check the h5py installation and the file then try again.\n'
            crash(msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def _read_types_lammps_hdf5(self):

        """
        get from lammps *.h5 file. the data should have been written using commands like:
            dump            POSITIONS all h5md ${dt_dump} pos.h5 position species
            dump_modify     POSITIONS sort id
        """

        self.types = np.zeros(self.num_atoms,dtype=int)

        try:
            import h5py
            with h5py.File(self.trajectory_file,'r') as in_db:

                # it is assumed that the types are consecutive integers starting at 1
                # so that 1 is subtracted to match the python index starting at 0 ...
                self.types[:] = in_db['particles/all/species/value'][0,:]-1

            # error check
            if self.num_atoms != self.types.size:
                msg = 'number of atoms in file dont match the input info.'
                msg += 'correct the configuration or check your trajectory file!\n'
                crash(msg)

        except Exception as _ex:
            msg = f'the hdf5 file\n  \'{self.trajectory_file}\'\ncould not be read!\n'
            msg += 'check the h5py installation and the file then try again.\n'
            crash(msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def _read_pos_lammps_hdf5(self,inds):

        """
        get from lammps *.h5 file. the data should have been written using commands like:
            dump            POSITIONS all h5md ${dt_dump} pos.h5 position species
            dump_modify     POSITIONS sort id
        """

        self.timers.start_timer('read_lammps_hdf5',units='s')

        try:
            import h5py
            with h5py.File(self.trajectory_file,'r') as in_db:

                if self.read_box:
                    self.box_vectors[0,0] = np.mean(
                            in_db['particles/all/box/edges/value'][inds,0],axis=0)
                    self.box_vectors[1,1] = np.mean(
                            in_db['particles/all/box/edges/value'][inds,1],axis=0)
                    self.box_vectors[2,2] = np.mean(
                            in_db['particles/all/box/edges/value'][inds,2],axis=0)

                # read positions
                self.pos[:,:,:] = in_db['particles/all/position/value'][inds,:,:]

        except Exception as _ex:
            msg = f'the hdf5 file\n  \'{self.trajectory_file}\'\ncould not be read!\n'
            msg += 'check the h5py installation and the file then try again.\n'
            crash(msg,_ex)

        self.timers.stop_timer('read_lammps_hdf5')

    # ----------------------------------------------------------------------------------------------

    def _read_pos_user_hdf5(self,inds):

        """
        get from user created hdf5 file
        """

        self.timers.start_timer('read_user_hdf5',units='s')

        if self.read_box:
            msg = '\nbox vectors must be given by user to use user_hdf5 positions (for now)'
            crash(msg)

        try:
            import h5py
            with h5py.File(self.trajectory_file,'r') as in_db:
                # read positions
                self.pos[:,:,:] = in_db['cartesian_pos'][inds,:,:]

        except Exception as _ex:
            msg = f'the hdf5 file\n  \'{self.trajectory_file}\'\ncould not be read!\n'
            msg += 'check the h5py installation and the file then try again.\n'
            crash(msg,_ex)

        self.timers.stop_timer('read_user_hdf5')

    # ----------------------------------------------------------------------------------------------

    def _read_pos_external(self,inds):

        """
        simply slice data from previously specified array
        """

        if self.read_box:
            msg = '\nbox vectors must be given by user to use external positions'
            crash(msg)

        self.pos[:,:,:] = self.external_pos[inds,:,:]

    # ----------------------------------------------------------------------------------------------


