#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                               !
#   ! this file is part of the 'pynamic-structure-factor' code      !
#   ! written by Ty Sterling at the University of Colorado          !
#   ! Boulder, advised by Dmitry Reznik.                            !
#   !                                                               !
#   ! the software calculates inelastic neutron dynamic structure   !
#   ! factors from molecular dynamics trajectories.                 !
#   !                                                               !
#   ! this is free software distrubuted under the GNU GPL v3 and    !
#   ! with no warrantee or garauntee of the results. you should     !
#   ! have recieved a copy of the new license with this software    !
#   ! if you do find bugs or have questions, dont hesitate to       !
#   ! write to the author at ty.sterling@colorado.edu               !
#   !                                                               !
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import os
import numpy as np
from mod_utils import PSF_exception, print_stdout


class input_variables:

    """
    object to hold all the input variables. sets defaults and parses input file to overwrite them
    if a keyword appears twice in the input file, only the 1st occurence is used. should change this
    in the future
    """ 

    # -----------------------------------------------------------------------------------------

    def __init__(self):
        """
        set the defaults
        """
        # this list is used to check that variables in input file are right. best practice is to add
        # any new key_words to this list
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
                          'compute_bragg',
                          'compute_sqw',
                          'compute_timeavg',
                          'parse_custom']

        self.doubles = [] # check if any keywords are duplicated

        self.traj_file      = 'pos.hdf5'    # *.hdf5 file holding MD positions
        self.output_dir     = 'sqw_output'  # where to write the *.hdf5 files holding sqw 
        self.outfile_prefix = 'sqw'         # prefix for output files, e.g. sqw_*_.hdf5
        self.save_progress  = False         # print output file for each process after each block

        # the default is to read the standard (using libhdf5) output form lammps.
        # the parse custom option is to allow to read different formats. modify the 
        # method if mod_io
        self.parse_custom   = False         # whether to read some other custom hdf5 heirarchy

        self.dt              = 1            # MD timestep in units of fs
        self.stride          = 32           # stride between when positions are printed
        self.total_steps     = 2**21        # total number of MD steps (NOT number of data points)
        self.num_atoms       = 4096         # number of atoms in the MD simulation
        self.supercell       = [[8,8,8]]    # size of supercell to (optionally) recalc. cell lengths
        self.lattice_vectors = [5.431,0,0,0,5.431,0,0,0,5.431] # only ortho cells for now

        self.unwrap_pos               = True    # unimpose minimum image convention
        self.recalculate_cell_lengths = True    # recalculate cell lens from box bounds in traj file
        self.b                        = [4.1491]    # neutron scatt. lengths
        self.Qpoints_file             = False       # read Q points from file (give name of file here)
        self.Qmin                     = [0,0,0]     # min of Q path to make 2d scan
        self.Qmax                     = [2,0,0]     # max of Q path to make 2d scan
        self.total_Qsteps             = 17          # number of steps in 2d scan

        # split the trajectory into 'blocks' and compute sqw for each. average the results of 
        # each block, sort of like an 'ensemble' average. can use unecessarily long trajectories
        # to improve statistics. the code automatically loops over the blocks. 
        self.num_blocks               = 1           # number of blocks to split

        # list of 'blocks'. e.g. if num_blocks = 4, this could be blocks = 0 1 2 3 or blocks = 0 etc.
        # defaults to all blocks, but you can pick 1 or 2 or ... useful for debugging. set num_blocks
        # to a large number and blocks = 0 so that the code only loops over 1 (short trajectory) block
        self.blocks                   = list(range(self.num_blocks))

        # compute 'bragg' scattering intensity and write it to '*_bragg.hdf5'
        # the bragg intensity is |<rho(Q)>|^2 (see Dove, appendix E.)
        # the total (time averaged) scattered intensity is <rho(Q)*rho(-Q)>
        # which is equal to the integral of S(Q,w) over all freqs. the bragg
        # intensity is defined as the total intensity minus the bragg intensity
        self.compute_bragg          = False
        # compute the timeaveraged total scattered intensity
        self.compute_timeavg        = False
        # compute the dynamical scattered intensity S(Q,w)
        self.compute_sqw            = True 

    # -----------------------------------------------------------------------------------------

    def parse_input(self,input_file):
        """
        read the input file 
        """
        # test that the input file exists/isnt broken
        try:
            with open(input_file,'r') as inp:
                self.input_txt = inp.readlines()
        except:
            message = f'input file \'{input_file}\' not found'
            raise PSF_exception(message)

        # check the key_words in the file, remove empty lines and comments
        self._check_file()

        # get the variables from file
        self.traj_file      = self._parse_str('traj_file',self.traj_file)   
        self.outfile_prefix = self._parse_str('outfile_prefix',self.outfile_prefix)
        self.output_dir     = self._parse_str('output_dir',self.output_dir)
        self.save_progress  = self._parse_bool('save_progress',self.save_progress)
        self.parse_custom   = self._parse_bool('parse_custom',self.parse_custom)
        self.dt          = self._parse_float('dt',self.dt)
        self.stride      = self._parse_int('stride',self.stride)
        self.total_steps = self._parse_int('total_steps',self.total_steps)
        self.num_atoms   = self._parse_int('num_atoms',self.num_atoms)
        self.supercell   = self._parse_int_list('supercell',self.supercell)
        self.lattice_vectors = self._parse_float_list('lattice_vectors',self.lattice_vectors)
        self.unwrap_pos               = self._parse_bool('unwrap_pos',self.unwrap_pos)
        self.recalculate_cell_lengths = self._parse_bool('recalculate_cell_lengths',
                                                        self.recalculate_cell_lengths)
        self.b         = self._parse_float_list('b',self.b)
        self.num_types = len(self.b)
        self.Qpoints_file = self._parse_str('Qpoints_file',self.Qpoints_file)
        self.Qmin         = self._parse_float_list('Qmin',self.Qmin)
        self.Qmax         = self._parse_float_list('Qmax',self.Qmax)
        self.total_Qsteps = self._parse_int('total_Qsteps',self.total_Qsteps)
        self.num_blocks   = self._parse_int('num_blocks',self.num_blocks)
        self.blocks       = list(range(self.num_blocks)) 
        self.blocks       = self._parse_int_list('blocks',self.blocks) 
        self.compute_bragg   = self._parse_bool('compute_bragg',self.compute_bragg)
        self.compute_timeavg = self._parse_bool('compute_timeavg',self.compute_timeavg)
        self.compute_sqw     = self._parse_bool('compute_sqw',self.compute_sqw)

        # check that the lattice vectors make sense
        try:
            self.lattice_vectors = np.array(self.lattice_vectors).reshape((3,3))
        except:
            message = 'lattice vectors seem wrong. should be a list of 9 floats with no commas'
            raise PSF_exception(message)
        # check that lattice vectors are ortho
        # the issue is that positions etc. are in cartesian coords with ortho boxes. different lattice
        # vectors should work, but i haven't tested it yet. it will be necessary to convert Q in 1/A 
        # to cartesian coordinates so that the vectorized multiplication done in mod_sqw._loop_over_blocks
        # works. 
        if (self.lattice_vectors[0,1] != 0 or self.lattice_vectors[0,2] != 0 or 
            self.lattice_vectors[1,0] != 0 or self.lattice_vectors[1,2] != 0  or
            self.lattice_vectors[2,0] != 0  or self.lattice_vectors[2,1] != 0): 
            message = 'only ortho. lattice vectors are currently supported. see comments in mod_invars'
            raise PSF_exception(message)

        # print the traj file
        message = f'reading trajectories from file \'{self.traj_file}\''
        print_stdout(message,msg_type='NOTE')

        # print the scattering lengths to file
        message = f' atom-type:  0    b: {self.b[0]: 2.4f}\n'
        for bb in range(1,self.num_types):
            message = message+f'  atom-type: {bb:2g}    b: {self.b[bb]: 2.4f}\n'
        print_stdout(message,msg_type='scattering lengths (b) in femtometers')

        # now convert b to Angstrom. 'intensities' are of order 1 this way. with FM, theyre 
        # of order 1e16. this is purely for convenience and doesnt change the physics, since the 
        # intensities are arbitrary units anyway
        for bb in range(self.num_types):
            self.b[bb] = self.b[bb]*1e-5

        # check that Q paths opts make sense
        if len(self.Qmin) != 3:
            message = f'variable Qmin should be a list of 3 floats'
            raise PSF_exception(message)
        if len(self.Qmax) != 3:
            message = f'variable Qmax should be a list of 3 floats'
            raise PSF_exception(message)

        # check that the requested blocks make sense
        if max(self.blocks) >= self.num_blocks or len(self.blocks) > self.num_blocks:
            message = f'variable blocks should be a list of the blocks to calculate'
            raise PSF_exception(message)

        # if the output dir. doesnt exist, create it
        if not os.path.exists(self.output_dir):
            message = f'creating directory \'{self.output_dir}\''
            print_stdout(message,msg_type='NOTE')
            os.mkdir(self.output_dir)

        # check that atleast one of compute_* is not False
        if not self.compute_sqw and not self.compute_timeavg and not self.compute_bragg:
            message = ('there is nothing to do! set atleast one of compute_sqw, \n compute_timeavg,'
                       ' or compute_bragg to 1 in the input file')
            raise PSF_exception(message)

        # check if traj file opens
        if not os.path.exists(self.traj_file):
            message = f'file \'{self.traj_file}\' not found'
            raise PSF_exception(message)

    # =======================================================================================
    # ------------------------------ private methods ----------------------------------------
    # =======================================================================================

    def _check_file(self):
        """
        check the key_words in input files and remove comments/blank lines
        """
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
                if key_word in self.doubles:
                    message = f'key word \'{key_word}\' appears more than once in the input file'
                    raise PSF_exception(message)
                self.doubles.append(key_word)
                input_txt.append(tmp_line)

    # -----------------------------------------------------------------------------------------

    def _parse_str(self,key_word,default):
        """
        get str varaible from file
        """
        return_value = default
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                return_value = str(return_value)
        return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_float(self,key_word,default):
        """
        get float varible from file
        """
        return_value = default
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
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
        """
        get int variable from file
        """
        return_value = default
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
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
        """
        get bool variable from file
        """
        return_value = default
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                try:
                    return_value = bool(int(return_value))
                except:
                    message = f'key word \'{key_word}\' seems wrongs.'
                    raise PSF_exception(message)
        return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_int_list(self,key_word,default):
        """
        get list of ints from file
        """
        return_value = default
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
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

    def _parse_float_list(self,key_word,default):
        """
        get list of floats from file
        """
        return_value = default
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                return_value = return_value.split()
                try:
                    return_value = [float(x) for x in return_value]
                except:
                    message = f'key word \'{key_word}\' seems wrongs.'
                    raise PSF_exception(message)
        return return_value

    # ------------------------------------------------------------------------------------------


