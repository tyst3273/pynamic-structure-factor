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

from psf.m_error import crash, check_file


class c_vasp_xdatcar_reader:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,input_files,number_of_atoms,constant_volume,block_size):

        """
        class to read number of steps, positions, etc from vasp5.x XDATCAR files
        """

        self.input_files = input_files
        self.num_files = len(input_files)
        self.number_of_atoms = number_of_atoms
        self.constant_volume = constant_volume
        self.block_size = block_size

        self._get_data_from_files()

    # ----------------------------------------------------------------------------------------------

    def _get_data_from_files(self):

        """
        get the number of steps in the files
        """

        self.lines_in_files = []
        self.types = []
        self.num_types = []
        self.num_each_type = []

        for _file in self.input_files:

            with open(_file,'r') as _if:

                # get number of lines in the file
                for _num_lines, _ in enumerate(_if):
                    pass
                self.lines_in_files.append(_num_lines)

                # go back to beginning of file
                _if.seek(0)

                # get types in file
                for ii in range(6):
                    _line = _if.readline()
                _types = _line.strip().split()
                self.types.append(_types)
                self.num_types.append(len(_types))
        
                # get number of each type
                _line = _if.readline().strip().split()
                _num_each_type = [int(_) for _ in _line]
                self.num_each_type.append(_num_each_type)

        # check that files are compatible
        if self.num_files > 1:
            _types = self.types[0]
            _num_types = self.num_types[0]
            _num_each = self.num_each_type[0]

            _msg = 'data in files are incompatible\n'
            for ii in range(1,self.num_files):
                if not _types == self.types[ii]:
                    crash(_msg)
                if not _num_types == self.num_types[ii]:
                    crash(_msg)
                if not _num_each == self.num_each_type[ii]:
                    crash(_msg)

        self.types = self.types[0]
        self.num_types = self.num_types[0]
        self.num_each_type = self.num_each_type[0]
        
        _num_atoms_in_file = sum(self.num_each_type)
        if not _num_atoms_in_file == self.number_of_atoms:
            _msg = 'number of atoms in files is wrong\n'
            crash(_msg)

    # ----------------------------------------------------------------------------------------------

    def get_steps_in_files(self):

        """
        get the number of steps in the files 
        """

        steps = []

        for ii in range(self.num_files):
            _n = self.lines_in_files[ii]

            # for NVT, lattice vectors only written to file at 1st step
            if self.constant_volume:
                _s = (_n-6)/(self.number_of_atoms+1)
                steps.append(int(_s))

            # for NPT, lattice vectors are written too
            else:
                _s = (_n+1)/(self.number_of_atoms+8)
                steps.append(int(_s))

        self.steps_in_files = steps

        return self.steps_in_files

    # ----------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------


class c_user_data:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,output_file='user_data.hdf5',input_files='XDATCAR',number_of_atoms=None,
                 file_format='vasp_xdatcar',constant_volume=True,block_size=None):

        """
        a module to read externally generated data from various codes (e.g. vasp, lammps, ...)
        in various formats. data are read from file in 'blocks' so that the whole file can be 
        parsed consecutively in chunks rather than reading the whole thing into memory all at once


        output_file is where the trajectory is written

        input_files is either a single file as a string or multiple files as a list of strings.
        if multiple files are given, they all have to have the same number of atoms 

        file_format specifies what methods to use to read the files. allowed values are 
            vasp_xdatcar  --  use this for VASP5.x XDATCAR files
            lammpstrj  -- usethis for lammpstrj trajectory files read by vmd
            ... add more if needed

        constant_volume tells wheter or not to read lattice vectors from file;
        FOR VASP: in NPT ensemble, simulation box fluctuates and supercell lattice
        vectors are writen to the file. set constant_volume=False for NPT.
        in NVT, they are not written so set constant_volume=True

        number_of_atoms is self explanatory... must be same for all files.

        block_size gives the maximum number of steps to be read/written at a time.
        the code will automatically divide up each file in the right way
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

        self.file_format = str(file_format)
        if self.file_format not in ['vasp_xdatcar','lammpstrj']:
            crash(f'unknown file_format type \'{self.file_format}\'\n')
        if self.file_format == 'vasp_xdatcar':
            self.reader = c_vasp_xdatcar_reader(self.input_files,self.number_of_atoms,
                                                self.constant_volume,self.block_size)

    # ----------------------------------------------------------------------------------------------

    def merge_files(self):

        """
        merge the files together into a single hdf5 database
        """
    
        # get number of steps in each file
        self.steps_in_files = self.reader.get_steps_in_files()
        print('steps in each file:\n',self.steps_in_files)
        self.total_steps = sum(self.steps_in_files)
        print('\ntotal number of steps:',self.total_steps)
       
        # get indices of steps to make blocks from
        self._get_block_inds()

        # get indices of lines in the file to read from/to
#        self.reader.get_file_inds(self.block_inds)

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
                self.block_inds.append(np.array([0,_steps],dtype=int))
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


# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    xdat = c_user_data(number_of_atoms=64,input_files=['npt/run_1/XDATCAR','npt/run_2/XDATCAR'],
            file_format='vasp_xdatcar',constant_volume=False,block_size=20)
    xdat.merge_files()

    xdat = c_user_data(number_of_atoms=64,input_files='nvt/XDATCAR',
            file_format='vasp_xdatcar',constant_volume=True)
    xdat.merge_files()











