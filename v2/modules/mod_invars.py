import os
import numpy as np
from mod_utils import PSF_exception, print_stdout


class input_variables:

    """
    object to hold all the input variables
    """ 

    # -----------------------------------------------------------------------------------------

    def __init__(self):
        """
        set the defaults
        """
        self.key_words = ['traj_file',
                          'output_dir',
                          'outfile_prefix',
                          'save_progress',
                          'dt',
                          'stride',
                          'total_steps',
                          'num_atoms',
                          'supercell',
                          'lattice_vectors',
                          'unwrap_pos',
                          'recalculate_cell_lengths',
                          'b',
                          'Qpoints_file',
                          'Qmin',
                          'Qmax',
                          'total_Qsteps',
                          'num_blocks',
                          'blocks',
                          'parse_lammps']

        self.traj_file      = 'pos.hdf5'
        self.output_dir     = 'sqw_output'
        self.outfile_prefix = 'sqw'
        self.save_progress  = False
        self.parse_lammps   = False

        self.dt              = 1e-15
        self.stride          = 32
        self.total_steps     = 2**21
        self.num_atoms       = 4096
        self.supercell       = [[8,8,8]]
        self.lattice_vectors = [5.431,0,0,0,5.431,0,0,0,5.431]

        self.unwrap_pos               = True
        self.recalculate_cell_lengths = True
        self.b                        = [4.1491e-5] 
        self.Qpoints_file             = False
        self.Qmin                     = [0,0,0]
        self.Qmax                     = [2,0,0]
        self.total_Qsteps             = 17
        self.num_blocks               = 1 
        self.blocks                   = list(range(self.num_blocks))

    # -----------------------------------------------------------------------------------------

    def parse_input(self,input_file):
        """
        read the input file 
        """
        try:
            with open(input_file,'r') as inp:
                self.input_txt = inp.readlines()
        except:
            message = f'input file \'{input_file}\' not found'
            raise PSF_exception(message)

        self._check_file()

        self.traj_file      = self._parse_str('traj_file',self.traj_file)

        self.output_dir     = self._parse_str('output_dir',self.output_dir)
        self.outfile_prefix = self._parse_str('outfile_prefix',self.outfile_prefix)
        if not os.path.exists(self.output_dir):
            message = f'creating directory \'{self.output_dir}\''
            print_stdout(message,msg_type='NOTE')
            os.mkdir(self.output_dir)

        self.save_progress  = self._parse_bool('save_progress',self.save_progress)
        self.parse_lammps   = self._parse_bool('parse_lammps',self.parse_lammps)

        self.dt          = self._parse_float('dt',self.dt)
        self.stride      = self._parse_int('stride',self.stride)
        self.total_steps = self._parse_int('total_steps',self.total_steps)
        self.num_atoms   = self._parse_int('num_atoms',self.num_atoms)
        self.supercell   = self._parse_int_vec('supercell',self.supercell)

        self.lattice_vectors = self._parse_float_vec('lattice_vectors',self.lattice_vectors)
        self.lattice_vectors = np.array(self.lattice_vectors).reshape((3,3))

        self.unwrap_pos               = self._parse_bool('unwrap_pos',self.unwrap_pos)
        self.recalculate_cell_lengths = self._parse_bool('recalculate_cell_lengths',
                                                        self.recalculate_cell_lengths)

        self.b         = self._parse_float_vec('b',self.b)
        self.num_types = len(self.b)
        message = ''
        for bb in range(self.num_types):
            this_b  = self.b[bb]
            message = message+f'  atom-type: {bb}    b: {this_b:2.4f}\n'
        print_stdout(message,msg_type='scattering lengths (b) in femtometers')

        # now convert b to Angstrom. 'intensities' are ~ 1 this way. with FM, theyre 
        # ~ 1e16
        for bb in range(self.num_types):
            self.b[bb] = self.b[bb]*1e-5

        self.Qpoints_file = self._parse_str('Qpoints_file',self.Qpoints_file)
        self.Qmin         = self._parse_float_vec('Qmin',self.Qmin)
        self.Qmax         = self._parse_float_vec('Qmax',self.Qmax)
        if len(self.Qmin) != 3:
            message = f'variable Qmin should be a list of 3 floats'
            raise PSF_exception(message)
        if len(self.Qmax) != 3:
            message = f'variable Qmax should be a list of 3 floats'
            raise PSF_exception(message)

        self.total_Qsteps = self._parse_int('total_Qsteps',self.total_Qsteps)
        self.num_blocks   = self._parse_int('num_blocks',self.num_blocks)
        self.blocks       = list(range(self.num_blocks)) 
        self.blocks       = self._parse_int_vec('blocks',self.blocks) 
        if max(self.blocks) >= self.num_blocks or len(self.blocks) > self.num_blocks:
            message = f'variable blocks should be a list of the blocks to calculate'
            raise PSF_exception(message)

    # =======================================================================================
    # ------------------------------ private methods ----------------------------------------
    # =======================================================================================

    def _check_file(self):
        input_txt = []
        for line in self.input_txt:
            if len(line.split()) == 0 or line.strip().startswith('#'):
                continue
            else:
                tmp_line = line.split('#')[0].strip()
                key_word = tmp_line.split('=')[0].strip()
                if key_word not in self.key_words:
                    message = f'key word \'{key_word}\' is unknown. check the input file'
                    raise PSF_exception(message)
                input_txt.append(tmp_line)

    # -----------------------------------------------------------------------------------------

    def _parse_str(self,key_word,default):
            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    return_value = str(return_value)
            return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_float(self,key_word,default):
            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    try:
                        return_value = float(return_value)
                    except:
                        message = f'key word \'{key_word}\' seems wrongs.'
                        raise PSF_exception(message)
            return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_int(self,key_word,default):
            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    try:
                        return_value = int(return_value)
                    except:
                        message = f'key word \'{key_word}\' seems wrongs.'
                        raise PSF_exception(message)
            return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_bool(self,key_word,default):
            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    try:
                        return_value = bool(int(return_value))
                    except:
                        message = f'key word \'{key_word}\' seems wrongs.'
                        raise PSF_exception(message)
            return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_int_vec(self,key_word,default):
            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    return_value = return_value.split()
                    try:
                        return_value = [int(x) for x in return_value]
                    except:
                        message = f'key word \'{key_word}\' seems wrongs.'
                        raise PSF_exception(message)
            return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_float_vec(self,key_word,default):
            return_value = default
            for line in self.input_txt:
                if line.strip().startswith(key_word):
                    return_value = line.split('=')[-1]
                    return_value = return_value.split('#')[0].strip()
                    return_value = return_value.split()
                    try:
                        return_value = [float(x) for x in return_value]
                    except:
                        message = f'key word \'{key_word}\' seems wrongs.'
                        raise PSF_exception(message)
            return return_value


