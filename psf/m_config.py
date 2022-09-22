
# system modules
import numpy as np
import importlib
import argparse
import os

# custom modules
import psf.m_import as m_import
from psf.m_error import check_file, crash


eps = 0.0001

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

        # get cmd line args
        description = 'command line args for \'PSF.py\''
        cmd_parser = argparse.ArgumentParser(description=description)

        # input file
        help_msg = 'input file for \'PSF.py\'. it is a python file that is imported as a module'
        cmd_parser.add_argument('-i','--input-file',default='input_params.py',help=help_msg)

        # get cmd line args
        cmd_args = cmd_parser.parse_args()

        # get input file from cmd line
        if input_file is None:
            self.input_file = cmd_args.input_file
        else:
            self.input_file = input_file

        # check if input file exists
        check_file(self.input_file) 

        # get defaults which are the allowed keywords... probably shouldt be hard-coded
        self.defaults = m_import.import_module('psf.defaults')

    # ----------------------------------------------------------------------------------------------

    def get_args_from_file(self):

        """
        get args from file
        """

        allowed = vars(self.defaults).keys()

        # load input file
        self.input = m_import.import_module_from_path(self.input_file)

        # check that all keywords are allowed
        for key in vars(self.input).keys():
            if key not in allowed:
                msg = f'the keyword \'{key}\' is unknown\n'
                crash(msg)

        print('\n*** input-args ***')
        print('echoing input args to the screen!\n')

        # now get the variables
        self._set_trajectory_format()
        self._set_trajectory_file()
        self._set_unwrap_trajectory()
        self._set_recalculate_box_lengths()
        self._set_calc_bragg()
        self._set_calc_diffuse()
        self._set_calc_sqw()
        self._set_lattice_vectors()
        self._set_atom_types()
        self._set_output_directory()
        self._set_output_prefix()
        self._set_md_time_step()
        self._set_md_num_steps()
        self._set_md_num_atoms()
        self._set_md_supercell_reps()
        self._set_num_trajectory_blocks()
        self._set_trajectory_blocks()
        self._set_experiment_type() 
        self._set_num_Qpoint_procs()
        self._set_Qpoints_option() 

        # could just read all of this stuff; this is just to make it lightweight 
        if self.Qpoints_option == 'mesh':
            self._set_Q_mesh_symmetry()
            self._set_Q_mesh_H()
            self._set_Q_mesh_K()
            self._set_Q_mesh_L()

            if self.Q_mesh_symmetry:
                self._set_basis_positions()

        elif self.Qpoints_option == 'text_file' \
                or self.Qpoints_option == 'mesh_file' \
                or self.Qpoints_option == 'write_mesh':
            self._set_Q_file() 
        elif self.Qpoints_option == 'path':
            self._set_Q_path_start()
            self._set_Q_path_end()
            self._set_Q_path_steps()

   
        # --- add new variables parsers here ---
        # ...
        # ...

    # ----------------------------------------------------------------------------------------------

    def _set_num_trajectory_blocks(self,num_trajectory_blocks=None):

        """
        get the attribute
        """

        if num_trajectory_blocks is None:
            self.num_trajectory_blocks = self._get_attr('num_trajectory_blocks')
        else:
            self.num_trajectory_blocks = num_trajectory_blocks
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

    def _set_trajectory_blocks(self,trajectory_blocks=None):

        """
        get the attribute
        """

        if trajectory_blocks is None:
            self.trajectory_blocks = self._get_attr('trajectory_blocks')
        else:
            self.trajectory_blocks = trajectory_blocks

        if self.trajectory_blocks is None:
            self.trajectory_blocks = np.arange(self.num_trajectory_blocks)
        else:
            msg = 'trajectory_blocks must be a list of ints'
            try:
                self.trajectory_blocks = np.array(self.trajectory_blocks,dtype=int)
            except:
                crash(msg)

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

    def _set_Q_path_steps(self,Q_path_steps=None):

        """
        get the attribute
        """

        if Q_path_steps is None:
            self.Q_path_steps = self._get_attr('Q_path_steps')
        else:
            self.Q_path_steps = Q_path_steps

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

    def _set_Q_path_end(self,Q_path_end=None):

        """
        get the attribute
        """

        if Q_path_end is None:
            self.Q_path_end = self._get_attr('Q_path_end')
        else:
            self.Q_path_end = Q_path_end

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

    def _set_Q_path_start(self,Q_path_start=None):

        """
        get the attribute
        """

        if Q_path_start is None:
            self.Q_path_start = self._get_attr('Q_path_start')
        else:
            self.Q_path_start = Q_path_start
        
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

    def _set_Q_file(self,Q_file=None):

        """
        get the attribute
        """

        if Q_file is None:
            self.Q_file = self._get_attr('Q_file')
        else:
            self.Q_file = Q_file
        self.Q_file = str(self.Q_file)
        self.Q_file = os.path.abspath(self.Q_file)

        # check if file exists
        check_file(self.Q_file)

        print(f'Q_file:\n  {self.Q_file}')

    # ----------------------------------------------------------------------------------------------

    def _set_Q_mesh_K(self,Q_mesh_K=None):

        """
        get the attribute
        """

        if Q_mesh_K is None:
            self.Q_mesh_K = self._get_attr('Q_mesh_K')
        else:
            self.Q_mesh_K = Q_mesh_K

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
            msg = f'Q_mesh_K:\n  {self.Q_mesh_K: 8.5f}'
        else:
            msg = 'Q_mesh_K:\n  '
            msg += f'{self.Q_mesh_K[0]: 8.5f} {self.Q_mesh_K[1]: 8.5f} {self.Q_mesh_K[2]:3g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_Q_mesh_L(self,Q_mesh_L=None):

        """
        get the attribute
        """

        if Q_mesh_L is None:
            self.Q_mesh_L = self._get_attr('Q_mesh_L')
        else:
            self.Q_mesh_L = Q_mesh_L

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
            msg = f'Q_mesh_L:\n  {self.Q_mesh_L: 8.5f}'
        else:
            msg = 'Q_mesh_L:\n  '
            msg += f'{self.Q_mesh_L[0]: 8.5f} {self.Q_mesh_L[1]: 8.5f} {self.Q_mesh_L[2]:3g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_Q_mesh_H(self,Q_mesh_H=None):

        """ 
        get the attribute 
        """

        if Q_mesh_H is None:
            self.Q_mesh_H = self._get_attr('Q_mesh_H')
        else:
            self.Q_mesh_H = Q_mesh_H

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
            msg = f'Q_mesh_H:\n  {self.Q_mesh_H: 8.5f}'
        else:
            msg = 'Q_mesh_H:\n  '
            msg += f'{self.Q_mesh_H[0]: 8.5f} {self.Q_mesh_H[1]: 8.5f} {self.Q_mesh_H[2]:3g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_Qpoints_option(self,Qpoints_option=None):

        """
        get the attribute
        """

        if Qpoints_option is None:
            self.Qpoints_option = self._get_attr('Qpoints_option')
        else:
            self.Qpoints_option = Qpoints_option
        self.Qpoints_option = str(self.Qpoints_option)

        msg = 'Qpoints_option seems wrong\n'
        if self.Qpoints_option not in ['text_file','mesh_file','write_mesh','mesh','path']:
            crash(msg)

        print(f'Qpoints_option:\n  {self.Qpoints_option}')

    # ----------------------------------------------------------------------------------------------

    def _set_experiment_type(self,experiment_type=None):

        """
        get the attribute
        """

        if experiment_type is None:
            self.experiment_type = self._get_attr('experiment_type')
        else:
            self.experiment_type = experiment_type
        self.experiment_type = str(self.experiment_type)

        msg = 'experiment_type seems wrong\n'
        if self.experiment_type not in ['neutrons','xrays']:
            crash(msg)

        print(f'experiment_type:\n  {self.experiment_type}')

    # ----------------------------------------------------------------------------------------------

    def _set_md_supercell_reps(self,md_supercell_reps=None):

        """
        get the attribute
        """

        if md_supercell_reps is None:
            self.md_supercell_reps = self._get_attr('md_supercell_reps')
        else:
            self.md_supercell_reps = md_supercell_reps

        msg = 'md_supercell_reps seems wrong\n'
        try:
            self.md_supercell_reps = np.array(self.md_supercell_reps,dtype=int)
        except:
            crash(msg)
        if self.md_supercell_reps.size != 3:
            crash(msg)

        if np.any(self.md_supercell_reps < 1):
            crash(msg)

        msg ='md_supercell_reps:\n  '
        for ii in range(3):
            msg = msg+f'{self.md_supercell_reps[ii]} '
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_md_num_atoms(self,md_num_atoms=None):

        """
        get the attribute
        """

        if md_num_atoms is None:
            self.md_num_atoms = self._get_attr('md_num_atoms')
        else:
            self.md_num_atoms = md_num_atoms
        self.md_num_atoms = int(self.md_num_atoms)

        msg = 'md_num_atoms must be > 0'
        if self.md_num_atoms <= 0:
            crash(msg)

        print(f'md_num_atoms:\n  {self.md_num_atoms}')

    # ----------------------------------------------------------------------------------------------

    def _set_md_num_steps(self,md_num_steps=None):

        """
        get the attribute
        """

        if md_num_steps is None:
            self.md_num_steps = self._get_attr('md_num_steps')
        else:
            self.md_num_steps = md_num_steps
        self.md_num_steps = int(self.md_num_steps)

        msg = 'md_num_steps must be > 0'
        if self.md_num_steps <= 0:
            crash(msg)

        print(f'md_num_steps:\n  {self.md_num_steps}')

    # ----------------------------------------------------------------------------------------------

    def _set_md_time_step(self,md_time_step=None):

        """
        get the attribute
        """

        if md_time_step is None:
            self.md_time_step = self._get_attr('md_time_step')
        else:
            self.md_time_step = md_time_step
        self.md_time_step = float(self.md_time_step)

        msg = 'md_time_step must be > 0'
        if self.md_time_step <= 0.0:
            crash(msg)

        print(f'md_time_step:\n  {self.md_time_step}')

    # ----------------------------------------------------------------------------------------------

    def _set_output_prefix(self,output_prefix=None):

        """
        get the attribute
        """

        if output_prefix is None:
            self.output_prefix = self._get_attr('output_prefix')
        else:
            self.output_prefix = output_prefix
        self.output_prefix = str(self.output_prefix+'_')
        print(f'output_prefix:\n  {self.output_prefix}')

    # ----------------------------------------------------------------------------------------------

    def _set_output_directory(self,output_directory=None):

        """
        get the attribute
        """

        if output_directory is None:
            self.output_directory = self._get_attr('output_directory')
        else:
            self.output_directory = output_directory

        # get the directory
        if self.output_directory is None:
            self.output_directory = os.getcwd()
        else:
            self.output_directory = os.path.abspath(self.output_directory)

        # make it if it doesnt exist
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)

        print(f'output_directory:\n  {self.output_directory}')

    # ----------------------------------------------------------------------------------------------

    def _set_calc_bragg(self,calc_bragg=None):

        """
        get the attribute
        """

        if calc_bragg is None:
            self.calc_bragg = self._get_attr('calc_bragg')
        else:
            self.calc_bragg = calc_bragg
        self.calc_bragg = bool(self.calc_bragg)
        print(f'calc_bragg:\n  {self.calc_bragg}')

    # ----------------------------------------------------------------------------------------------

    def _set_calc_diffuse(self,calc_diffuse=None):

        """
        get the attribute
        """

        if calc_diffuse is None:
            self.calc_diffuse = self._get_attr('calc_diffuse')
        else:
            self.calc_diffuse = calc_diffuse
        self.calc_diffuse = bool(self.calc_diffuse)
        print(f'calc_diffuse:\n  {self.calc_diffuse}')

    # ----------------------------------------------------------------------------------------------

    def _set_calc_sqw(self,calc_sqw=None):

        """
        get the attribute
        """

        if calc_sqw is None:
            self.calc_sqw = self._get_attr('calc_sqw')
        else:
            self.calc_sqw = calc_sqw
        self.calc_sqw = bool(self.calc_sqw)
        print(f'calc_sqw:\n  {self.calc_sqw}')

    # ----------------------------------------------------------------------------------------------

    def _set_recalculate_box_lengths(self,recalculate_box_lengths=None):

        """
        get the attribute
        """

        if recalculate_box_lengths is None:
            self.recalculate_box_lengths = self._get_attr('recalculate_box_lengths')
        else:
            self.recalculate_box_lengths = recalculate_box_lengths
        self.recalculate_box_lengths = bool(self.recalculate_box_lengths)

        # --- devloper stuff ---
        if self.recalculate_box_lengths:
            msg = '\'recalculate_box_lengths\' is disabled right now. if you need this\n' \
                  'functionality, contact the author at --- ty.sterling@colorado.edu ---\n'
            crash(msg)
        # ----------------------

        #print(f'recalculate_box_lengths:\n  {self.recalculate_box_lengths}')

    # ----------------------------------------------------------------------------------------------

    def _set_unwrap_trajectory(self,unwrap_trajectory=None):

        """
        get the attribute
        """

        if unwrap_trajectory is None:
            self.unwrap_trajectory = self._get_attr('unwrap_trajectory')
        else:
            self.unwrap_trajectory = unwrap_trajectory
        self.unwrap_trajectory = bool(self.unwrap_trajectory)
        print(f'unwrap_trajectory:\n  {self.unwrap_trajectory}')

    # ----------------------------------------------------------------------------------------------

    def _set_trajectory_file(self,trajectory_file=None):

        """
        get the attribute
        """

        if trajectory_file is None:
            self.trajectory_file = self._get_attr('trajectory_file')
        else:
            self.trajectory_file = trajectory_file
        self.trajectory_file = str(self.trajectory_file)
        self.trajectory_file = os.path.abspath(self.trajectory_file)

        # check if file exists
        check_file(self.trajectory_file)

        print(f'trajectory_file:\n  {self.trajectory_file}')

    # ----------------------------------------------------------------------------------------------

    def _set_trajectory_format(self,trajectory_format=None):

        """
        get the attribute
        """

        if trajectory_format is None:
            self.trajectory_format = self._get_attr('trajectory_format')
        else:
            self.trajectory_format = trajectory_format
        self.trajectory_format = str(self.trajectory_format)
        
        msg = 'trajectory_format \'{self.trajectory_format}\' is unknown'
        if self.trajectory_format not in ['lammps_hdf5']:
            crash(msg)

        print(f'trajectory_format:\n  {self.trajectory_format}')

    # ----------------------------------------------------------------------------------------------

    def _set_lattice_vectors(self,lattice_vectors=None):
        
        """
        get the attribute
        """
    
        if lattice_vectors is None:
            self.lattice_vectors = self._get_attr('lattice_vectors')
        else:
            self.lattice_vectors = lattice_vectors

        # check the shape
        msg = f'\'lattice_vectors\' seems wrong\n'
        try:
            self.lattice_vectors = np.array(self.lattice_vectors,dtype=float)
        except:
            crash(msg)
        if self.lattice_vectors.size != 9:
            crash(msg)
        else:
            self.lattice_vectors.shape = [3,3]

        # check if lattice vectors are numerically orthorhombic
        self.ortho_lattice_vectors = True
        for ii in range(3):
            for jj in range(3):
                if ii == jj:
                    continue
                else:
                    if np.round(self.lattice_vectors[ii,jj],6) != 0.0:
                        self.ortho_lattice_vectors = False

        if self.ortho_lattice_vectors:
            msg = f'orthorhombic_lattice_vectors:\n  {self.ortho_lattice_vectors}'
            print(msg)
        else:
            msg = 'only orthogonal (i.e. diagonal) lattice_vectors allowed for now\n'
            msg += f'if this functionality is really needed, contact the author'
            msg += ' at\n --- ty.sterling@colorado.edu\n'
            crash(msg)

        # echo to the info file
        msg = f'lattice_vectors: (Angstrom)'
        for ii in range(3):
            msg += '\n'
            for jj in range(3):
                msg += f'{self.lattice_vectors[ii,jj]: 10.6f} '
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_basis_positions(self,basis_positions=None):

        """
        get the attribute
        """

        if basis_positions is None:
            self.basis_positions = self._get_attr('basis_positions')
        else:
            self.basis_positions = basis_positions

        # if None, won't use them to reduce symmetry
        if self.basis_positions is None:
            self.use_basis_positions = False
            msg = 'basis_positions:\n  None (not using symmetry)'
            print(msg)

        else:
            # check the shape
            msg = f'\'basis_positions\' seems wrong\n'
            try:
                self.basis_positions = np.array(self.basis_positions,dtype=float)
            except:
                crash(msg)
            if self.basis_positions.shape[1] != 3:
                crash(msg)
            if self.basis_positions.shape[0] != self.num_basis_atoms:
                crash(msg)

            # echo to the info file
            msg = f'basis_positions: (crystal coords.)'
            for ii in range(self.num_basis_atoms):
                msg = msg+'\n'
                for jj in range(3):
                    msg = msg+f'{self.basis_positions[ii,jj]: 10.6f} '
            print(msg)
    
    # ----------------------------------------------------------------------------------------------

    def _set_atom_types(self,atom_types=None):

        """
        get the attribute

        IMPORTANT: note how the types are parsed; if the same type is repeated after other 
        are defined, e.g. ['Si','Si','C','Si','Ge'], the result for unique_types will be 
        'Si', 'C', 'Ge'.  the unique types (i.e. those that 'atom_files' are mapped to) 
        will be used in the order *new* types are found in 'atom_types' arg
        """
        
        if atom_types is None:
            self.atom_types = self._get_attr('atom_types')
        else:
            self.atom_types = atom_types

        # number of atoms in unitcell
        self.num_basis_atoms = len(self.atom_types)

        # find unique types; can use np.unique because dont want to sort it
        self.unique_types = []
        for atom in self.atom_types:
            if atom in self.unique_types:
                continue
            else:
                self.unique_types.append(atom)
        self.num_types = len(self.unique_types)

        # get type map
        self.type_map = []
        for ii in range(self.num_basis_atoms):
            self.type_map.append(self.unique_types.index(self.atom_types[ii]))

        # echo to the info file
        msg = 'atom_types:\n  '
        for ii in range(self.num_basis_atoms):
            msg = msg+f'{self.atom_types[ii]} '
        print(msg)
        msg = f'num_basis_atoms:\n  {self.num_basis_atoms:g}\n'
        msg = msg+f'num_types:\n  {self.num_types:g}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_Q_mesh_symmetry(self,Q_mesh_symmetry=None):

        """
        get the attribute
        """

        if Q_mesh_symmetry is None:
            self.Q_mesh_symmetry = self._get_attr('Q_mesh_symmetry')
        else:
            self.Q_mesh_symmetry = Q_mesh_symmetry

        # check if allowed
        try:
            self.Q_mesh_symmetry = bool(self.Q_mesh_symmetry)
        except:
            msg = '\'Q_mesh_symmetry\' seems wrong'
            crash(msg)

        # echo to the info file
        msg = f'Q_mesh_symmetry:\n  {self.Q_mesh_symmetry}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _set_num_Qpoint_procs(self,num_Qpoint_procs=None):

        """
        get the attribute
        """

        if num_Qpoint_procs is None:
            self.num_Qpoint_procs = self._get_attr('num_Qpoint_procs')
        else:
            self.num_Qpoint_procs = num_Qpoint_proc1s
        
        # check if allowed
        msg = f'\'num_Qpoint_procs\' must be an integer that is > 0\n'
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
        check if input file has the attr we want. get it if so. 
        otherwise use default
        """

        # check if input file
        if hasattr(self.input,attr):
            return getattr(self.input,attr)

        # otherwise use default
        else:
            return getattr(self.defaults,attr)

    # ----------------------------------------------------------------------------------------------



