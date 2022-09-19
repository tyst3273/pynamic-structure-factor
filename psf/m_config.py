
# system modules
import numpy as np
import importlib
import argparse
import os

# custom modules
import psf.m_import as m_import
from psf.m_error import check_file, crash


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
        print(allowed)

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
        self._set_verbosity()

        self._set_lattice_vectors()
        self._set_atom_types()
        self._set_basis_positions()

        self._set_trajectory_format()
        self._set_trajectory_file()
        self._set_unwrap_trajectory()
        self._set_recalculate_box_lengths()

        self._set_calc_bragg()
        self._set_calc_timeavg()
        self._set_calc_sqw()
        self._set_output_directory('./out/')
        self._set_output_prefix()

        self._set_md_time_step()
        self._set_md_num_steps()
        self._set_md_num_atoms()
        self._set_md_supercell_reps()

        self._set_experiment_type() 
        self._set_Qpoints_option() 
        self._set_Q_mesh_symmetry()
        self._set_Q_mesh() 
#        self._set_Q_mesh_H()
#        self._set_Q_mesh_K()
#        self._set_Q_mesh_L()
#        self._set_Q_file() 
#        self._set_Q_path_start()
#        self._set_Q_path_end()
#        self._set_Q_path_steps()
   
        self._set_num_Qpoint_procs()

        # --- add new variables parsers here ---
        # ...
        # ...

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
        if self.Qpoints_option not in ['file','path','mesh']:
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

    def _set_calc_timeavg(self,calc_timeavg=None):

        """
        get the attribute
        """

        if calc_timeavg is None:
            self.calc_timeavg = self._get_attr('calc_timeavg')
        else:
            self.calc_timeavg = calc_timeavg
        self.calc_timeavg = bool(self.calc_timeavg)
        print(f'calc_timeavg:\n  {self.calc_timeavg}')

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
        print(f'recalculate_box_lengths:\n  {self.recalculate_box_lengths}')

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

        # echo to the info file
        msg = f'lattice_vectors: (Angstrom)'
        for ii in range(3):
            msg = msg+'\n'
            for jj in range(3):
                msg = msg+f'{self.lattice_vectors[ii,jj]: 10.6f} '
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

    def _set_Q_mesh(self,Q_mesh=None):

        """
        get the attribute
        """

        if Q_mesh is None:
            self.Q_mesh = self._get_attr('Q_mesh')
        else:
            self.Q_mesh = Q_mesh

        # check the shape
        msg = f'\'Q_mesh\' should be 3 integers > 0\n'
        try:
            self.Q_mesh = np.array(self.Q_mesh,dtype=int)
        except:
            crash(msg)
        if self.Q_mesh.size != 3:
            crash(msg)
        if np.any(self.Q_mesh < 1):
            crash(msg)

        # echo to the info file
        msg = f'Q_mesh:\n  {self.Q_mesh[0]:g} {self.Q_mesh[1]:g} {self.Q_mesh[2]:g}'
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

    def _set_verbosity(self,verbosity=None):

        """
        get the attribute
        """

        if verbosity is None:
            self.verbosity = self._get_attr('verbosity')
        else:
            self.verbosity = verbosity

        # check if allowed
        try:
            self.verbosity = int(self.verbosity)
        except:
            msg = f'\'verbosity\' must 0, 1, or 2'
            crash(msg)
        if self.verbosity not in [0,1,2]:
            msg = f'\'verbosity\' must 0, 1, or 2'
            crash(msg)

        # echo to the info file
        msg = f'verbosity:\n  {self.verbosity:g}'
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



