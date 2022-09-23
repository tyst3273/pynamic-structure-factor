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
        template
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config
        self.timers = timers

        # the trajectory file
        self.num_atoms = self.config.md_num_atoms
        self.trajectory_format = self.config.trajectory_format
        self.trajectory_file = self.config.trajectory_file
        self.md_time_step = self.config.md_time_step
    
        self.unwrap_trajectory = self.config.unwrap_trajectory

        # set up the 'blocks' of indices for calculating on
        self._get_block_inds()
    
        # the times for time correlation functions
        self.time = np.arange(0,self.num_block_steps,self.md_time_step)

        # get atom types once and for all
        self._read_types_lammps_hdf5()

        # WRAPPED positions in cartesian coords in Angstrom
        self.pos = np.zeros((self.num_block_steps,self.num_atoms,3))

     # ----------------------------------------------------------------------------------------------

    def _get_block_inds(self):

        """
        setup indices for block averaging
        """

        _md_steps = self.config.md_num_steps
        _num_blocks = self.config.num_trajectory_blocks

        # averging is done in frequency space; having blocks of data that are on
        # different frequency grids will require interpolation to average together
        # i want to avoid doing that so i crop the dataset so that all blocks
        # have same frequency
        if _md_steps % _num_blocks != 0:
            msg = '\n*** warning! ***\n'
            msg += 'md_num_steps must be divisible by num_trajectory_blocks\n'
            msg += 'discarding the extra timesteps\n'
            print(msg)

        self.blocks = self.config.trajectory_blocks
        
        _rem = _md_steps % _num_blocks
        self.num_block_avg = self.config.num_block_avg
        self.num_steps_used = _md_steps-_rem
        self.num_block_steps = self.num_steps_used // _num_blocks

        # lo/hi inds to index the blocks from the trajectory file 
        # note that full trajectory is indexed and only the blocks that are 
        # requested are used
        _blocks = np.arange(_num_blocks)
        self.block_inds = np.zeros((_num_blocks,2),dtype=int)
        self.block_inds[:,0] = _blocks*self.num_block_steps
        self.block_inds[:,1] = (_blocks+1)*self.num_block_steps

        msg = f'\n*** blocks ***\nnum_block_avg: {self.num_block_avg}\n'
        msg += f'num_block_steps: {self.num_block_steps}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def get_atom_types(self):

        """
        get an array that holds types of atom. only do for 1st step ...

        it is MANDATORY that data are sorted in the exact same way each timestep... the
        code assumes the order of the types of atoms doesnt change over time, so only calculates
        the scattering info at the first timestep
        """

        if self.trajectory_format == 'lammps_hdf5':
            self._read_types_lammps_hdf5(_inds)

    # ----------------------------------------------------------------------------------------------

    def read_trajectory_block(self,block_index=0):

        """
        get a block of trajectory data from the file and preprocess it as much as possible
        """

        _inds = self.block_inds[block_index]

        if self.trajectory_format == 'lammps_hdf5':
            self._read_pos_lammps_hdf5(_inds)            

        if self.unwrap_trajectory:
            self._unwrap_positions()

     # ----------------------------------------------------------------------------------------------

    def _unwrap_positions(self):

        """
        un-apply minimum image convention so that positions dont jump discontinuously
        by a full box length. this wont be memory effecient but should be fast... ish

        see 'note about unwrapping' in the README
        """

        self.timers.start_timer('unwrap_positions',units='s')

        msg = '\nunwrapping positions!\n'
        print(msg)

        # get box lengths... assumes the simulation cell and underlying unitcell are orthorhombic
        lx = self.comm.lattice.lattice_vectors[0,0]*self.config.md_supercell_reps[0]
        ly = self.comm.lattice.lattice_vectors[1,1]*self.config.md_supercell_reps[1]
        lz = self.comm.lattice.lattice_vectors[2,2]*self.config.md_supercell_reps[2]

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
                self.types[:] = in_db['particles']['all']['species']['value'][0,:]-1

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

    def _read_pos_lammps_hdf5(self,_inds):

        """
        get from lammps *.h5 file. the data should have been written using commands like:
            dump            POSITIONS all h5md ${dt_dump} pos.h5 position species
            dump_modify     POSITIONS sort id
        """

        self.timers.start_timer('read_lammps_hdf5',units='s')

        try:

            import h5py
            with h5py.File(self.trajectory_file,'r') as in_db:

                # not currently used ...
                """
                self.box_lengths[0] = np.mean(
                    in_db['particles']['all']['box']['edges']['value'][inds[0]:inds[1],0],axis=0)
                self.box_lengths[1] = np.mean(
                    in_db['particles']['all']['box']['edges']['value'][inds[0]:inds[1],1],axis=0)
                self.box_lengths[2] = np.mean(
                    in_db['particles']['all']['box']['edges']['value'][inds[0]:inds[1],2],axis=0)
                """
                # read positions
                self.pos[:,:,:] = in_db['particles']['all']['position']['value'] \
                                    [_inds[0]:_inds[1],:,:]

        except Exception as _ex:
            msg = f'the hdf5 file\n  \'{self.trajectory_file}\'\ncould not be read!\n'
            msg += 'check the h5py installation and the file then try again.\n'
            crash(msg,_ex)

        self.timers.stop_timer('read_lammps_hdf5')

    # ----------------------------------------------------------------------------------------------


