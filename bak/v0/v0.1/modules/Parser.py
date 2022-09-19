import os
import numpy as np
import h5py 
from FileIO import print_stdout
from ParalUtils import ParalExcept

class parser:

    def __init__(self):

        """
        set the defaults for all options. 
        """
        
        # list of keywords to check that none are wrong
        self.key_words = ['pos_dir','traj_file','output_dir','outfile_prefix','save_progress',
                'dt','stride','total_steps','num_atoms','supercell','lattice_vectors','unwrap_pos',
                'recalculate_cell_lengths','b','Qpoints_file','Qmin','Qmax','total_Qsteps',
                'num_blocks','blocks','parse_lammps']

        # ============= defaults =================

        # I/O
        self.pos_dir = os.getcwd()
        self.traj_file = os.path.join(self.pos_dir,'pos.hdf5')
        self.output_dir = 'sqe_output'
        self.outfile_prefix = 'sqe'
        self.save_progress = False
        self.parse_lammps = False

        # MD params
        self.dt = 1e-15
        self.stride = 32
        self.total_steps = 2**21
        self.num_atoms = 4096
        self.supercell = [[8,8,8]]
        self.lattice_vectors = [5.431,0,0,0,5.431,0,0,0,5.431]

        # calc options
        self.unwrap_pos = True
        self.recalculate_cell_lengths = True
        self.b = [4.1491e-5]
        self.Qpoints_file = False
        self.Qmin = [0,0,0]
        self.Qmax = [2,0,0]
        self.total_Qsteps = 17
        self.num_blocks = 1
        self.blocks = list(range(self.num_blocks)) 



    def parse(self,input_file):

        """
        parse the input file for all keywords. if the option isnt in the input file, return 
        the default.
        """

        try:
            with open(input_file,'r') as inp:
                self.input_txt = inp.readlines()
        except:
            message = f'input file \'{input_file}\' not found'
            raise ParalExcept(message)

        # ========== check keywords in file =========
        self._check_file()

        # ================ I/O ======================
        self.pos_dir = self._parse_str('pos_dir',self.pos_dir)
        self.traj_file = self._parse_str('traj_file',self.traj_file)
        self.traj_file = os.path.join(self.pos_dir,self.traj_file)
        self.output_dir = self._parse_str('output_dir',self.output_dir)
        self.outfile_prefix = self._parse_str('outfile_prefix',self.outfile_prefix)
        self.save_progress = self._parse_bool('save_progress',self.save_progress)
        self.parse_lammps = self._parse_bool('parse_lammps',self.parse_lammps)

        # =============== MD params =====================
        self.dt = self._parse_float('dt',self.dt)
        self.stride = self._parse_int('stride',self.stride)
        self.total_steps = self._parse_int('total_steps',self.total_steps)
        self.num_atoms = self._parse_int('num_atoms',self.num_atoms)
        self.supercell = self._parse_int_vec('supercell',self.supercell)

        # this is a hack to keep format consistent with the rest of my code
        self.lattice_vectors = self._parse_float_vec('lattice_vectors',self.lattice_vectors)
        lattice_vectors = [[0,0,0],[0,0,0],[0,0,0]]
        ind = 0
        for ii in range(3):
            for jj in range(3):
                lattice_vectors[ii][jj] = self.lattice_vectors[ind]
                ind = ind+1
        self.lattice_vectors = lattice_vectors

        # =================== calc options ==================
        self.unwrap_pos = self._parse_bool('unwrap_pos',self.unwrap_pos)
        self.recalculate_cell_lengths = self._parse_bool('recalculate_cell_lengths',
                                                        self.recalculate_cell_lengths)

        # this is to keep format consistent with the rest of the code
        self.b = self._parse_float_vec('b',self.b)
        self.num_types = len(self.b)
        b = {} # i dunno why i made this a dict originally, but here it is. 
        for bb in range(self.num_types):
            b[f'{bb+1}'] = self.b[bb]
        self.b = b

        self.Qpoints_file = self._parse_str('Qpoints_file',self.Qpoints_file)
        self.Qmin = self._parse_float_vec('Qmin',self.Qmin)
        self.Qmax = self._parse_float_vec('Qmax',self.Qmax)
        self.total_Qsteps = self._parse_int('total_Qsteps',self.total_Qsteps)
        self.num_blocks = self._parse_int('num_blocks',self.num_blocks)
        self.blocks = list(range(self.num_blocks)) # update default
        self.blocks = self._parse_int_vec('blocks',self.blocks) # overwrite default IF in the file.


    # ====================================================
    # ---------------- private methods -------------------
    # ====================================================
    
    def _check_file(self):

        for line in self.input_txt:

            # skip blank lines and comment lines
            if len(line.split()) == 0 or line.strip().startswith('#'):
                continue

            else:
                key_word = line.strip().split('=')[0].strip()
                if key_word not in self.key_words:
                    message = f'\'{key_word}\' is not one of the expected variables. check the input file'
                    raise ParalExcept(message=message)

    def _parse_str(self,key_word,default):

            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    return_value = str(return_value)

            return return_value

    def _parse_float(self,key_word,default):

            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    return_value = float(return_value)

            return return_value

    def _parse_int(self,key_word,default):

            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    return_value = int(return_value)

            return return_value

    def _parse_bool(self,key_word,default):

            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    return_value = bool(int(return_value))

            return return_value

    def _parse_int_vec(self,key_word,default):

            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    return_value = return_value.split()
                    return_value = [int(x) for x in return_value]
                    return_value = return_value

            return return_value

    def _parse_float_vec(self,key_word,default):

            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    return_value = return_value.split()
                    return_value = [float(x) for x in return_value]
                    return_value = return_value

            return return_value






