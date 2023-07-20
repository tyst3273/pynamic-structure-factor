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

# system modules
import numpy as np
import argparse
import os

# custom modules
import psf.m_import as m_import
from psf.m_error import check_file, crash
from psf.m_empty import c_empty


# --------------------------------------------------------------------------------------------------

def get_input_files():

    """
    get input file(s) from command line args. default is None (dont use).
    """

    # get cmd line args
    description = 'command line args for \'PSF.py\''
    cmd_parser = argparse.ArgumentParser(description=description)

    # input file
    help_msg = 'input file for \'PSF.py\'. it is a python file that is imported as a module\n'
    help_msg += 'to batch run input files, give multiple in the order you want to run them'
    cmd_parser.add_argument('-i','--input-files',default=[None],help=help_msg,nargs='*')

    # get cmd line args
    cmd_args = cmd_parser.parse_args()
    input_files = cmd_args.input_files

    return input_files

# --------------------------------------------------------------------------------------------------


class c_config:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,input_file=None):

        """
        get cmd line args and check for input file. 

        NOTE: none of the attribute in this class should be overwritten. 
        it is ONLY for getting input args at execution time. everything else should be 
        kept in more appropriate places. nothing should be assigned here either. 
        it is "read-only"
        """

        # get defaults which are the allowed keywords... probably shouldt be hard-coded
        self.defaults = m_import.import_module('psf.defaults')
        self.allowed = vars(self.defaults).keys()

        print(f'\n*** input-args ***')

        # get input file from cmd line
        if input_file is not None:
            check_file(input_file)
            self.input = m_import.import_module_from_path(input_file)
            print(f'reading input from {input_file}\n')
        else:
            self.input = c_empty()
            print(f'using defaults for non-keyword args\n')

        # check that all keywords in input file are allowed
        for key in vars(self.input).keys():
            if key not in self.allowed:
                msg = f'the keyword \'{key}\' is unknown\n'
                crash(msg)

    # ----------------------------------------------------------------------------------------------

    def set_config(self,**kwargs):

        """
        get args from file; if args are given in function call, they are put in self.kwargs and 
        overwrite whatever is read from input and default.  priority is arg > input > default.
        """

        # these are args passed into api
        self.kwargs = kwargs

        # check that all keywords args are allowed
        for key in self.kwargs.keys():
            if key not in self.allowed:
                msg = f'the keyword \'{key}\' is unknown\n'
                crash(msg)

        # now get the variables
        self._set_trajectory_format()
        if not self.trajectory_format in ['external']:
            self._set_trajectory_file()

        self._set_unwrap_trajectory()
        self._set_calc_sqw()

        self._set_lattice_vectors()

        if self.unwrap_trajectory:
            self._set_box_vectors()

        self._set_atom_types()
        self._set_output_directory()
        self._set_output_prefix()
        self._set_md_time_step()
        self._set_md_num_steps()
        self._set_md_num_atoms()
        self._set_trajectory_stride()
        self._set_trajectory_skip()
        self._set_trajectory_trim()
        self._set_num_trajectory_blocks()
        self._set_trajectory_blocks()
        self._set_experiment_type() 
        self._set_num_Qpoint_procs()
        self._set_Qpoints_option() 

        # could just read all of this stuff; this is just to make it lightweight 
        if self.Qpoints_option == 'mesh':
            self._set_Q_mesh_H()
            self._set_Q_mesh_K()
            self._set_Q_mesh_L()

        elif self.Qpoints_option == 'file':
            self._set_Q_file() 

        elif self.Qpoints_option == 'path':
            self._set_Q_path_start()
            self._set_Q_path_end()
            self._set_Q_path_steps()
            
        # --- add new variables parsers here ---
        # ...
        # ...

    # ----------------------------------------------------------------------------------------------

    def _set_trajectory_stride(self):

        """
        get the attribute
        """

        self._get_attr('trajectory_stride')

        try:
            self.trajectory_stride = int(self.trajectory_stride)
        except:
            msg = 'trajectory_stride must be an int'
            crash(msg)

        msg = 'trajectory_stride must be >= 1'
        if self.trajectory_stride < 1:
            crash(msg)

        print(f'trajectory_stride:\n  {self.trajectory_stride}')

    # ----------------------------------------------------------------------------------------------

    def _set_trajectory_skip(self):

        """
        get the attribute
        """

        self._get_attr('trajectory_skip')

        try:
            self.trajectory_skip = int(self.trajectory_skip)
        except:
            msg = 'trajectory_skip must be an int'
            crash(msg)

        msg = 'trajectory_skip must be >= 0'
        if self.trajectory_skip < 0:
            crash(msg)

        print(f'trajectory_skip:\n  {self.trajectory_skip}')    

    # ----------------------------------------------------------------------------------------------

    def _set_trajectory_trim(self):

        """
        get the attribute
        """

        self._get_attr('trajectory_trim')

        try:
            self.trajectory_trim = int(self.trajectory_trim)
        except:
            msg = 'trajectory_trim must be an int'
            crash(msg)

        msg = 'trajectory_trim must be >= 0'
        if self.trajectory_trim < 0:
            crash(msg)

        print(f'trajectory_trim:\n  {self.trajectory_trim}')

    # ----------------------------------------------------------------------------------------------

    def _set_num_trajectory_blocks(self):

        """
        get the attribute
        """

        self._get_attr('num_trajectory_blocks')

        try:
            self.num_trajectory_blocks = int(self.num_trajectory_blocks)
        except:
            msg = 'num_trajectory_blocks must be an int'
            crash(msg)

        msg = 'num_trajectory_blocks must be >= 1'
        if self.num_trajectory_blocks < 1:
            crash(msg)
    
        print(f'num_trajectory_blocks:\n  {self.num_trajectory_blocks}')

    # ----------------------------------------------------------------------------------------------

    def _set_trajectory_blocks(self):

        """
        get the attribute
        """

        self._get_attr('trajectory_blocks')

        if self.trajectory_blocks is None:
            self.trajectory_blocks = np.arange(self.num_trajectory_blocks)
        else:
            msg = 'trajectory_blocks must be a list of ints'
            try:
                self.trajectory_blocks = np.array(self.trajectory_blocks,dtype=int)
            except:
                crash(msg)

        if self.trajectory_blocks.ndim == 0:
            self.trajectory_blocks.shape = [1,]
        if len(self.trajectory_blocks.shape) != 1:
            crash(msg)

        # check that it makes sense
        self.num_block_avg = self.trajectory_blocks.size
        if self.trajectory_blocks.min() < 0 \
            or self.trajectory_blocks.max() >= self.num_trajectory_blocks:
            crash(msg)

        msg = 'trajectory_blocks:\n  '
        for ii in range(self.num_block_avg):
            msg += f'{self.trajectory_blocks[ii]:g} '
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_Q_path_steps(self):

        """
        get the attribute
        """

        self._get_attr('Q_path_steps')

        msg = 'Q_path_steps seems wrong\n'
        try:
            self.Q_path_steps = np.array(self.Q_path_steps,dtype=int)
        except:
            crash(msg)

        if len(self.Q_path_steps.shape) != 1:
            crash(msg)

        # error check
        msg = 'Q_path_steps, Q_path_start, and Q_path_end all should have same number of paths\n'
        self.num_Q_paths = self.Q_path_steps.size
        if self.Q_path_start.shape[0] != self.num_Q_paths:
            crash(msg)
        if self.Q_path_end.shape[0] != self.num_Q_paths:
            crash(msg)
        
        # print to screen
        msg = 'Q_path:' 
        for ii in range(self.num_Q_paths):
            msg = msg+f'\n  {self.Q_path_start[ii,0]:8.5f} {self.Q_path_start[ii,1]: 8.5f}' \
                      f' {self.Q_path_start[ii,2]:8.5f}  ==>  '
            msg = msg+f'{self.Q_path_end[ii,0]:8.5f} {self.Q_path_end[ii,1]:8.5f}' \
                      f' {self.Q_path_end[ii,2]:8.5f}    {self.Q_path_steps[ii]:3g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_Q_path_end(self):

        """
        get the attribute
        """

        self._get_attr('Q_path_end')

        msg = 'Q_path_end seems wrong\n'
        try:
            self.Q_path_end = np.array(self.Q_path_end,dtype=float)
        except:
            crash(msg)

        if self.Q_path_end.size == 3:
            self.Q_path_end.shape = [1,3]

        if self.Q_path_end.shape[1] != 3:
            crash(msg)

        # info is printed in _set_Q_path_steps()

    # ----------------------------------------------------------------------------------------------

    def _set_Q_path_start(self):

        """
        get the attribute
        """

        self._get_attr('Q_path_start')
        
        msg = 'Q_path_start seems wrong\n'
        try:
            self.Q_path_start = np.array(self.Q_path_start,dtype=float)
        except:
            crash(msg)

        if self.Q_path_start.size == 3:
            self.Q_path_start.shape = [1,3]

        if self.Q_path_start.shape[1] != 3:
            crash(msg)

        # info is printed in _set_Q_path_steps()

    # ----------------------------------------------------------------------------------------------

    def _set_Q_file(self):

        """
        get the attribute
        """

        self._get_attr('Q_file')
            
        self.Q_file = str(self.Q_file)
        self.Q_file = os.path.abspath(self.Q_file)

        # check if file exists
        check_file(self.Q_file)

        print(f'Q_file:\n  {self.Q_file}')

    # ----------------------------------------------------------------------------------------------

    def _set_Q_mesh_K(self):

        """
        get the attribute
        """

        self._get_attr('Q_mesh_K')

        msg = 'Q_mesh_K seems wrong\n'
        try:
            self.Q_mesh_K = np.array(self.Q_mesh_K,dtype=float)
        except:
            crash(msg)
        if not self.Q_mesh_K.size in [1,3]:
            crash(msg)
        if self.Q_mesh_K.size == 3:
            if self.Q_mesh_K[2] < 1:
                crash(msg)
        # if only 1 num, print it
        if self.Q_mesh_K.size == 1:
            if self.Q_mesh_K.ndim == 1:
                self.Q_mesh_K = self.Q_mesh_K[0]
            msg = f'Q_mesh_K:\n  {self.Q_mesh_K: 8.5f}'
        else:
            msg = 'Q_mesh_K:\n  '
            msg += f'{self.Q_mesh_K[0]: 8.5f} {self.Q_mesh_K[1]: 8.5f} {self.Q_mesh_K[2]:3g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_Q_mesh_L(self):

        """
        get the attribute
        """

        self._get_attr('Q_mesh_L')

        msg = 'Q_mesh_L seems wrong\n'
        try:
            self.Q_mesh_L = np.array(self.Q_mesh_L,dtype=float)
        except:
            crash(msg)
        if not self.Q_mesh_L.size in [1,3]:
            crash(msg)
        if self.Q_mesh_L.size == 3:
            if self.Q_mesh_L[2] < 1:
                crash(msg)

        # if only 1 num, print it
        if self.Q_mesh_L.size == 1:
            if self.Q_mesh_L.ndim == 1:
                self.Q_mesh_L = self.Q_mesh_L[0]
            msg = f'Q_mesh_L:\n  {self.Q_mesh_L: 8.5f}'
        else:
            msg = 'Q_mesh_L:\n  '
            msg += f'{self.Q_mesh_L[0]: 8.5f} {self.Q_mesh_L[1]: 8.5f} {self.Q_mesh_L[2]:3g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_Q_mesh_H(self):

        """ 
        get the attribute 
        """

        self._get_attr('Q_mesh_H')

        msg = 'Q_mesh_H seems wrong\n'
        try:
            self.Q_mesh_H = np.array(self.Q_mesh_H,dtype=float)
        except:
            crash(msg)
        if not self.Q_mesh_H.size in [1,3]:
            crash(msg)
        if self.Q_mesh_H.size == 3:
            if self.Q_mesh_H[2] < 1:
                crash(msg)

        # if only 1 num, print it
        if self.Q_mesh_H.size == 1:
            if self.Q_mesh_H.ndim == 1:
                self.Q_mesh_H = self.Q_mesh_H[0]
            msg = f'Q_mesh_H:\n  {self.Q_mesh_H: 8.5f}'
        else:
            msg = 'Q_mesh_H:\n  '
            msg += f'{self.Q_mesh_H[0]: 8.5f} {self.Q_mesh_H[1]: 8.5f} {self.Q_mesh_H[2]:3g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_Qpoints_option(self):

        """
        get the attribute
        """

        self._get_attr('Qpoints_option')

        self.Qpoints_option = str(self.Qpoints_option)

        msg = 'Qpoints_option seems wrong\n'
        if self.Qpoints_option not in ['file','mesh','path']:
            crash(msg)

        print(f'Qpoints_option:\n  {self.Qpoints_option}')

    # ----------------------------------------------------------------------------------------------

    def _set_experiment_type(self):

        """
        get the attribute
        """

        self._get_attr('experiment_type')

        self.experiment_type = str(self.experiment_type)

        msg = 'experiment_type seems wrong\n'
        if self.experiment_type not in ['neutrons','xrays']:
            crash(msg)

        print(f'experiment_type:\n  {self.experiment_type}')

    # ----------------------------------------------------------------------------------------------

    def _set_md_num_atoms(self):

        """
        get the attribute
        """

        self._get_attr('md_num_atoms')

        self.md_num_atoms = int(self.md_num_atoms)

        msg = 'md_num_atoms must be > 0'
        if self.md_num_atoms <= 0:
            crash(msg)

        print(f'md_num_atoms:\n  {self.md_num_atoms}')

    # ----------------------------------------------------------------------------------------------

    def _set_md_num_steps(self):

        """
        get the attribute
        """

        self._get_attr('md_num_steps')

        self.md_num_steps = int(self.md_num_steps)

        msg = 'md_num_steps must be > 0'
        if self.md_num_steps <= 0:
            crash(msg)

        print(f'md_num_steps:\n  {self.md_num_steps}')

    # ----------------------------------------------------------------------------------------------

    def _set_md_time_step(self):

        """
        get the attribute
        """

        self._get_attr('md_time_step')
            
        self.md_time_step = float(self.md_time_step)

        msg = 'md_time_step must be > 0'
        if self.md_time_step <= 0.0:
            crash(msg)

        print(f'md_time_step:\n  {self.md_time_step}')

    # ----------------------------------------------------------------------------------------------

    def _set_output_prefix(self):

        """
        get the attribute
        """

        self._get_attr('output_prefix')

        if self.output_prefix is None:
            return 

        self.output_prefix = str(self.output_prefix)
        print(f'output_prefix:\n  {self.output_prefix}')

    # ----------------------------------------------------------------------------------------------

    def _set_output_directory(self):

        """
        get the attribute
        """

        self._get_attr('output_directory')

        # get the directory
        if self.output_directory is None:
            self.output_directory = os.getcwd()
        else:
            self.output_directory = os.path.abspath(self.output_directory)

        print(f'output_directory:\n  {self.output_directory}')

    # ----------------------------------------------------------------------------------------------

    def _set_calc_sqw(self):

        """
        get the attribute
        """

        self._get_attr('calc_sqw')
            
        self.calc_sqw = bool(self.calc_sqw)
        print(f'calc_sqw:\n  {self.calc_sqw}')

    # ----------------------------------------------------------------------------------------------

    def _set_unwrap_trajectory(self):

        """
        get the attribute
        """

        self._get_attr('unwrap_trajectory')

        self.unwrap_trajectory = bool(self.unwrap_trajectory)
        print(f'unwrap_trajectory:\n  {self.unwrap_trajectory}')

    # ----------------------------------------------------------------------------------------------

    def _set_trajectory_file(self):

        """
        get the attribute
        """

        self._get_attr('trajectory_file')
            
        self.trajectory_file = str(self.trajectory_file)
        self.trajectory_file = os.path.abspath(self.trajectory_file)

        # check if file exists
        check_file(self.trajectory_file)

        print(f'trajectory_file:\n  {self.trajectory_file}')

    # ----------------------------------------------------------------------------------------------

    def _set_trajectory_format(self):

        """
        get the attribute
        """

        self._get_attr('trajectory_format')

        msg = 'trajectory_format \'{self.trajectory_format}\' is unknown'
        if self.trajectory_format not in ['lammps_hdf5','external','user_hdf5']:
            crash(msg)

        print(f'trajectory_format:\n  {self.trajectory_format}')

    # ----------------------------------------------------------------------------------------------

    def _set_lattice_vectors(self):
        
        """
        get the attribute
        """
    
        self._get_attr('lattice_vectors')

        # check the shape
        msg = '\'lattice_vectors\' seems wrong\n'
        try:
            self.lattice_vectors = np.array(self.lattice_vectors,dtype=float)
        except:
            crash(msg)
        if self.lattice_vectors.size != 9:
            crash(msg)
        else:
            self.lattice_vectors.shape = [3,3]

        # echo to the info file
        msg = 'lattice_vectors: (Angstrom)'
        for ii in range(3):
            msg += '\n'
            for jj in range(3):
                msg += f'{self.lattice_vectors[ii,jj]: 10.6f} '
        print(msg)

    # ----------------------------------------------------------------------------------------------
    
    def _set_box_vectors(self):

        """
        get the attribute
        """

        self._get_attr('box_vectors')

        if self.box_vectors is None:
            msg = 'box_vectors:\n read from file\n'
        return

        # check the shape
        msg = '\'box_vectors\' seems wrong\n'
        try:
            self.box_vectors = np.array(self.box_vectors,dtype=float)
        except:
            crash(msg)
        if self.box_vectors.size != 9:
            crash(msg)
        else:
            self.box_vectors.shape = [3,3]

        # echo to the info file
        msg = 'box_vectors: (Angstrom)'
        for ii in range(3):
            msg += '\n'
            for jj in range(3):
                msg += f'{self.box_vectors[ii,jj]: 10.6f} '
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_atom_types(self):

        """
        get the attribute.

        one string per type in the trajectory data. same string can be repeated for multiple 
        types
        """
        
        self._get_attr('atom_types')

        self.num_types = len(self.atom_types)

        # echo to the info file
        msg = 'atom_types:\n  '
        for ii in range(len(self.atom_types)):
            msg += f'{self.atom_types[ii]} '
        msg += f'\nnum_types:\n  {self.num_types:g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_num_Qpoint_procs(self):

        """
        get the attribute
        """

        self._get_attr('num_Qpoint_procs')
        
        # check if allowed
        msg = '\'num_Qpoint_procs\' must be an integer that is > 0\n'
        try:
            self.num_Qpoint_procs = int(self.num_Qpoint_procs)
        except:
            crash(msg)
        if self.num_Qpoint_procs <= 0:
            crash(msg)

        # echo to the info file
        msg = f'num_Qpoint_procs:\n  {self.num_Qpoint_procs:g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _get_attr(self,attr):

        """
        check if given as arg in kwargs. if so, use arg. if not, check if given in
        input file. if not, use default.
        """

        # check if arg is in kwargs
        if attr in self.kwargs.keys():
            arg = self.kwargs[attr]

        # otherwise get from input file or defaults
        else:

            # check if in input file
            if hasattr(self.input,attr):
                arg = getattr(self.input,attr)

            # otherwise use default
            else:
                arg = getattr(self.defaults,attr)

        setattr(self,attr,arg)

    # ----------------------------------------------------------------------------------------------

