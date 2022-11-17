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
import h5py
from timeit import default_timer

from psf.m_error import crash, check_file



class c_lammpstrj_reader:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,input_file,number_of_atoms,reduced,block_size):

        """
        class to read number of steps, positions, etc from lammpstrj files
        note that constant_volume has no meaning for this class
        """

        self.input_file = input_file
        self.number_of_atoms = number_of_atoms
        self.reduced = reduced
        self.block_size = block_size

        self._get_metadata()

    # ----------------------------------------------------------------------------------------------

    def _get_metadata(self):

        """
        get the number of steps etc. from the file
        """

        with open(self.input_file,'r') as _if:

            # get number of lines in the file
            for _num_lines, _ in enumerate(_if):
                pass
            _num_lines += 1

            # go back to beginning of file
            _if.seek(0)

            # get types in file
            for ii in range(9):
                _if.readline()

            _types = []
            for ii in range(self.number_of_atoms):
                _ = _if.readline().strip().split()
                _types.append(_[1])
            _unique = []
            self.types = []
            for ii in range(self.number_of_atoms):
                _t = _types[ii]
                if _t not in _unique:
                    _unique.append(_t)
                _i = _unique.index(_t)
                self.types.append(_i)

        # types start from 1
        self.types = [x+1 for x in self.types]           

        self.num_types = len(self.types)

        self.steps_in_file = _num_lines/(self.number_of_atoms+9)
        self.steps_in_file = int(self.steps_in_file)

    # ----------------------------------------------------------------------------------------------

    def open_file(self):

        """
        keep file open so we dont have to find the location we left each time we read from it
        """

        self.f = open(self.input_file,'r')

    # ----------------------------------------------------------------------------------------------

    def close_file(self):

        """
        close the file when done
        """

        self.f.close()

    # ----------------------------------------------------------------------------------------------

    def get_trajectory(self,num_steps,step_0=False):

        """
        get the sets of positions and box vectors from the file. num_steps is the number of steps
        to read (they are read consecutively since file is kept open)
        """

        self.pos = np.zeros((num_steps,self.number_of_atoms,3),dtype=float)
        self.box = np.zeros((num_steps,3,3),dtype=float)

        self._get_trj(num_steps)

        if not self.reduced:
            self.pos = self._get_reduced_coords(self.pos,self.box)

        return self.pos, self.box

    # ----------------------------------------------------------------------------------------------

    def _get_reduced_coords(self,pos,box_vecs):

        """
        go from cartesian coords to reduced coords
        """

        red = np.zeros(pos.shape)
        _num_steps = red.shape[0]
        for ii in range(_num_steps):
            _vecs = np.linalg.inv(box_vecs[ii,...])
            for jj in range(self.number_of_atoms):
                red[ii,jj,:] = pos[ii,jj,0]*_vecs[0,:]+ \
                               pos[ii,jj,1]*_vecs[1,:]+ \
                               pos[ii,jj,2]*_vecs[2,:]

        return red

    # ----------------------------------------------------------------------------------------------

    def _get_trj(self,num_steps):

        """
        get positions in reduced coords from file. since npt, read box vectors every step
        """

        for ii in range(num_steps):

            for jj in range(5):
                self.f.readline()

            _ = [float(x) for x in self.f.readline().strip().split()]
            self.box[ii,0,0] = _[1]-_[0]
            _ = [float(x) for x in self.f.readline().strip().split()]
            self.box[ii,1,1] = _[1]-_[0]
            _ = [float(x) for x in self.f.readline().strip().split()]
            self.box[ii,2,2] = _[1]-_[0]

            self.f.readline()

            # get positions
            for jj in range(self.number_of_atoms):
                self.pos[ii,jj,:] = self.f.readline().strip().split()[2:]

    # ----------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------


class c_vasp_xdatcar_reader:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,input_file,number_of_atoms,constant_volume,block_size):

        """
        class to read number of steps, positions, etc from vasp5.x XDATCAR files
        """

        self.input_file = input_file
        self.number_of_atoms = number_of_atoms
        self.constant_volume = constant_volume
        self.block_size = block_size

        self._get_metadata()

    # ----------------------------------------------------------------------------------------------

    def _get_metadata(self):

        """
        get the number of steps etc. from the file
        """

        with open(self.input_file,'r') as _if:

            # get number of lines in the file
            for _num_lines, _ in enumerate(_if):
                pass
            _num_lines += 1

            # go back to beginning of file
            _if.seek(0)

            # get types in file
            for ii in range(6):
                _line = _if.readline()
            _types = _line.strip().split()

            self.num_types = len(_types)
        
            # get number of each type
            _line = _if.readline().strip().split()
            _num_each_type = [int(_) for _ in _line]
            self.num_each_type = _num_each_type

        # for NVT, lattice vectors only written to file at 1st step
        if self.constant_volume:
            self.steps_in_file = (_num_lines-7)/(self.number_of_atoms+1)
        # for NPT, lattice vectors are written too
        else:
            self.steps_in_file = _num_lines/(self.number_of_atoms+8)
        self.steps_in_file = int(self.steps_in_file)
 
        self.types = []
        for ii in range(self.num_types):
            self.types.extend([ii]*self.num_each_type[ii])
        self.types = [x+1 for x in self.types]
    
    # ----------------------------------------------------------------------------------------------

    def open_file(self):

        """
        keep file open so we dont have to find the location we left each time we read from it
        """
        
        self.f = open(self.input_file,'r')

    # ----------------------------------------------------------------------------------------------

    def close_file(self):

        """
        close the file when done
        """

        self.f.close()

    # ----------------------------------------------------------------------------------------------

    def get_trajectory(self,num_steps,step_0=False):

        """
        get the sets of positions and box vectors from the file. num_steps is the number of steps
        to read (they are read consecutively since file is kept open)
        """

        self.pos = np.zeros((num_steps,self.number_of_atoms,3),dtype=float)
        self.box = np.zeros((num_steps,3,3),dtype=float)

        if self.constant_volume:
            self._get_nvt(num_steps,step_0)
        else:
            self._get_npt(num_steps)

        return self.pos, self.box

    # ----------------------------------------------------------------------------------------------

    def _get_npt(self,num_steps):

        """
        get positions in reduced coords from file. since npt, read box vectors every step
        """

        for ii in range(num_steps):

            self.f.readline()
            _s = float(self.f.readline().strip())
            self.box[ii,0,:] = self.f.readline().strip().split()
            self.box[ii,1,:] = self.f.readline().strip().split()
            self.box[ii,2,:] = self.f.readline().strip().split()
            self.box[ii,:,:] *= _s
            self.f.readline()
            self.f.readline()
            self.f.readline()

            # get positions
            for jj in range(self.number_of_atoms):
                self.pos[ii,jj,:] = self.f.readline().strip().split()

    # ----------------------------------------------------------------------------------------------

    def _get_nvt(self,num_steps,step_0):

        """
        get positions in reduced coords from file. since nvt, return the 1st step box vectors 
        instead of reading them. if 1st step, read the box vectors
        """

        # get the box vectors part if nvt and 1st step
        if step_0:
            self.box_0 = np.zeros((3,3))
            self.f.readline()
            _s = float(self.f.readline().strip())
            self.box_0[0,:] = self.f.readline().strip().split()
            self.box_0[1,:] = self.f.readline().strip().split()
            self.box_0[2,:] = self.f.readline().strip().split()
            self.box_0 *= _s
            self.f.readline()
            self.f.readline()

        for ii in range(num_steps):

            # skip 'Direct configuration = *' comment line
            self.f.readline()

            # assign box
            self.box[ii,:,:] = self.box_0[:,:]

            # get positions
            for jj in range(self.number_of_atoms):
                self.pos[ii,jj,:] = self.f.readline().strip().split()

    # ----------------------------------------------------------------------------------------------






# --------------------------------------------------------------------------------------------------

class c_user_data:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,output_file='user_data.hdf5',input_files='XDATCAR',number_of_atoms=None,
                 file_format='vasp_xdatcar',constant_volume=True,block_size=None,reduced=False):

        """
        a module to read externally generated data from various codes (e.g. vasp, lammps, ...)
        in various formats. data are read from file in 'blocks' so that the whole file can be 
        parsed consecutively in chunks rather than reading the whole thing into memory all at once

        output_file is where the trajectory is written

        input_files is either a single file as a string or multiple files as a list of strings.
        if multiple files are given, they all have to have the same number of atoms with atoms 
        written to the file in exactly the same order

        file_format specifies what methods to use to read the files. allowed values are 
            vasp_xdatcar  --  use this for VASP5.x XDATCAR files
            lammpstrj  -- use this for lammpstrj trajectory files that can be read by vmd
            ... add more if needed

        constant_volume tells wheter or not to read lattice vectors from file;
        FOR VASP: in NPT ensemble, simulation box fluctuates and supercell lattice
        vectors are writen to the file. set constant_volume=False for NPT.
        in NVT, they are not written so set constant_volume=True

        number_of_atoms is self explanatory... must be same for all files.

        block_size gives the maximum number of steps to be read/written at a time.
        the code will automatically divide up each file in the right way

        reduced tells wheter or not to expect data in crystal coords. set to False for
        cartesian coords

        """

        self.output_file = str(output_file)
        print(f'\nwriting to file:\n  \'{self.output_file}\'')

        if not isinstance(input_files,list): 
            input_files = [input_files]
        self.input_files = [str(_) for _ in input_files]            
        print('\nreading from files:\n',self.input_files)
        self.num_files = len(input_files)

        _msg = 'block_size must be an integer greater than 0\n'
        if block_size is not None:
            self.block_size = int(block_size)
            if self.block_size < 1:
                crash(_msg)
        else:
            self.block_size = block_size
            print('\nreading each file all at once!')

        _msg = 'number_of_atoms must be an integer greater than 0\n'
        if number_of_atoms is None:
            crash(_msg)
        self.number_of_atoms = int(number_of_atoms)
        if self.number_of_atoms < 1:
            crash(_msg)

        self.constant_volume = bool(constant_volume)
        self.reduced = bool(reduced)

        self.file_format = str(file_format)
        if self.file_format not in ['vasp_xdatcar','lammpstrj']:
            crash(f'unknown file_format type \'{self.file_format}\'\n')

        # set up list of 'readers' that get data from the file on command
        if self.file_format == 'vasp_xdatcar':
            self.readers = [c_vasp_xdatcar_reader(_f,self.number_of_atoms,
                            self.constant_volume,self.block_size) for _f in self.input_files]
        if self.file_format == 'lammpstrj':
            self.readers = [c_lammpstrj_reader(_f,self.number_of_atoms,
                            self.reduced,self.block_size) for _f in self.input_files]

    # ----------------------------------------------------------------------------------------------

    def merge_files(self):

        """
        merge the files together into a single hdf5 database
        """
    
        # get number of steps in each file
        self.steps_in_files = []
        for _r in self.readers:
            self.steps_in_files.append(_r.steps_in_file)
        print('\nsteps in each file:\n',self.steps_in_files)

        self.total_steps = sum(self.steps_in_files)
        print('\ntotal number of steps:',self.total_steps)
       
        # get indices of steps to make blocks from
        self._get_block_inds()

        # set up the hdf5 file
        with h5py.File(self.output_file,'w') as _odb:

            _odb.create_dataset('types',shape=(self.number_of_atoms),dtype=int)
            _odb.create_dataset('cartesian_pos',shape=(self.total_steps,self.number_of_atoms,3)
                    ,dtype=float)
            _odb.create_dataset('box_vectors',shape=(self.total_steps,3,3),dtype=float)

            # get first step from file for unwrapping the trajectory
            _r = self.readers[0]
            _r.open_file()
            pos_0, _ = _r.get_trajectory(num_steps=1,step_0=True)  
            _odb['types'][:] = _r.types[:]
            _r.close_file()

            _file_ind = 0
            for ii, _r in enumerate(self.readers):

                print(f'\nnow on file {ii}')
                _r.open_file()
                
                # loop over blocks and get the data from the files
                _block_inds = self.block_inds[ii]
                print(f'there are {_block_inds.shape[0]} blocks\n')
                for bb in range(_block_inds.shape[0]):  

                    _start_time = default_timer()
                    print(f'now on block {bb}')

                    if _block_inds[bb,0] == 0:
                        step_0 = True
                    else:
                        step_0 = False

                    # get steps from file
                    num_steps = _block_inds[bb,1]-_block_inds[bb,0]
                    pos, box = _r.get_trajectory(num_steps,step_0)

                    # unwrap and go to cartesian coords
                    pos = self._unwrap_pos(pos,pos_0)
                    pos = self._get_cartesian_pos(pos,box)

                    # writ this block to the file
                    _odb['cartesian_pos'][_file_ind:_file_ind+num_steps,...] = pos[...]
                    _odb['box_vectors'][_file_ind:_file_ind+num_steps,...] = box[...]

                    _file_ind += num_steps

                    _end_time = default_timer()
                    _elapsed = (_end_time-_start_time)
                    print(f'elapsed time for block: {_elapsed:10.6f} [s]\n')

                # close the text file
                _r.close_file()

    # ----------------------------------------------------------------------------------------------

    def _get_cartesian_pos(self,pos,box_vecs):

        """
        go from reduced to cartesian coords
        """

        cart = np.zeros(pos.shape)
        _num_steps = cart.shape[0]
        for ii in range(_num_steps):
            for jj in range(self.number_of_atoms):
                cart[ii,jj,:] = pos[ii,jj,0]*box_vecs[ii,0,:]+ \
                                pos[ii,jj,1]*box_vecs[ii,1,:]+ \
                                pos[ii,jj,2]*box_vecs[ii,2,:]

        return cart

    # ----------------------------------------------------------------------------------------------

    def _unwrap_pos(self,pos,pos_0):

        """
        unwrap positions back to the same convention as in step 0
        """

        num_steps = pos.shape[0]

        # build an array to be added to the positions to shift them 
        shift = pos-np.tile(pos_0[0,:,:].reshape((1,self.number_of_atoms,3)),
                reps=[num_steps,1,1])

        # check whether to shift atoms along 1st box vector
        dr = -(shift[:,:,0] > 1/2).astype(int)
        dr = dr+(shift[:,:,0] <= -1/2).astype(int)
        shift[:,:,0] = np.copy(dr)

        # check whether to shift atoms in the y direction  
        dr = -(shift[:,:,1] > 1/2).astype(int)
        dr = dr+(shift[:,:,1] <= -1/2).astype(int)
        shift[:,:,1] = np.copy(dr)

        # check whether to shift atoms in the z direction   
        dr = -(shift[:,:,2] > 1/2).astype(int)
        dr = dr+(shift[:,:,2] <= -1/2).astype(int)
        shift[:,:,2] = np.copy(dr)

        # apply the shift    
        return pos+shift

    # ----------------------------------------------------------------------------------------------

    def _get_block_inds(self):

        """
        divide up the steps into 'blocks' that are at most block_size steps long. 
        """

        _block_size = self.block_size

        self.block_inds = []

        for ii in range(self.num_files):
            _steps = self.steps_in_files[ii]
            
            if _block_size is None or _block_size >= _steps:
                self.block_inds.append(np.array([[0,_steps]],dtype=int))
                continue

            _n = int(np.ceil(_steps/_block_size))
            _r = int(_steps % _block_size)

            _b = np.arange(_n)
            _blocks = np.zeros((_n,2),dtype=int)
            _blocks[:,0] = _b*_block_size
            _blocks[:,1] = (_b+1)*_block_size

            if _r != 0:
                _blocks[-1,1] = _blocks[-1,1]-(_block_size-_r)
            self.block_inds.append(_blocks)

    # ----------------------------------------------------------------------------------------------

    def write_lammpstrj(self,trj_file='pos.lammpstrj',steps_to_write=None):
        
        """
        read the data from the hdf5 file into a lammpstrj text file. do NOT do this if the 
        trajectory file is huge...

        ITEM: TIMESTEP
            %d (timestep number)
        ITEM: NUMBER OF ATOMS
            %d (number of atoms)
        ITEM: BOX BOUNDS pp pp pp
            %f %f (xlo, xhi)
            %f %f (ylo, yhi)
            %f %f (zlo, zhi)
        ITEM: ATOMS id type x y z
            %d %d %f %f %f  (atomid, type, x-, y-, z, coordinates)

        """

        with h5py.File(self.output_file,'r') as _db, open(trj_file,'w') as _fo:
            
            _types = _db['types'][...]
            _shape = _db['cartesian_pos'].shape
            _num_steps = _shape[0]
            _num_atoms = _shape[1]

            if not steps_to_write is None:
                _num_steps = int(steps_to_write)
            
            for ii in range(_num_steps):

                _box = _db['box_vectors'][ii,...]
                _l0 = _box[0,0]; _l1 = _box[1,1]; _l2 = _box[2,2] 

                _ = 0.0
                _fo.write(f'ITEM: TIMESTEP\n{ii:g}\n')
                _fo.write(f'ITEM: NUMBER OF ATOMS\n{_num_atoms:g}\n')
                _fo.write('ITEM: BOX BOUNDS pp pp pp\n')
                _fo.write(f'{_:10.6f}  {_l0:10.6f}\n')
                _fo.write(f'{_:10.6f}  {_l1:10.6f}\n')
                _fo.write(f'{_:10.6f}  {_l2:10.6f}\n')
                _fo.write('ITEM: ATOMS id type x y z\n')
                
                _pos = _db['cartesian_pos'][ii,...]
                for jj in range(_num_atoms):
                    _fo.write(f' {jj+1:6g} {_types[jj]:3g} {_pos[jj,0]:12.9f}' \
                              f' {_pos[jj,1]:12.9f} {_pos[jj,2]:12.9f}\n')

    # ----------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    """
    this stuff is for debugging but is also an example of how to use the classes above
    """

    xdat = c_user_data(number_of_atoms=64,output_file='npt.hdf5',
            input_files=['vasp_sample/npt/run_1/XDATCAR','vasp_sample/npt/run_2/XDATCAR'],
            file_format='vasp_xdatcar',constant_volume=False,block_size=20)
    xdat.merge_files()
    xdat.write_lammpstrj('vasp_sample/npt.lammpstrj')

    xdat = c_user_data(number_of_atoms=64,output_file='nvt.hdf5',
            input_files='vasp_sample/nvt/XDATCAR',file_format='vasp_xdatcar',constant_volume=True)
    xdat.merge_files()
    xdat.write_lammpstrj('vasp_sample/nvt.lammpstrj')

    trj = c_user_data(number_of_atoms=64,output_file='npt_lmp.hdf5',
            input_files=['vasp_sample/npt.lammpstrj'],
            file_format='lammpstrj',block_size=20,reduced=False)
    trj.merge_files()

    trj = c_user_data(number_of_atoms=64,output_file='nvt_lmp.hdf5',
            input_files='vasp_sample/nvt.lammpstrj',file_format='lammpstrj',reduced=False)
    trj.merge_files()








