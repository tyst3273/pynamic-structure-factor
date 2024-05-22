#!/home/ty/anaconda3/bin/python

"""
Tool to read lammps text output and convert to format that can be read by my code 
i.e. trajectory_format='user_hdf5'.

I assume that the trajectory was written using the commands
    dump          pos all custom ${dump_freq} pos.dat id type xu yu zu
    dump_modify   pos sort id

It is **essential** that you used the unwrapped positions (xu, yu, zu) and sort by id!

This isn't parallelized (since parallel reads from the text file are apparently not faster).
It simply reads one step at a time, writing each to the hdf5 file. It isn't fast, but it is 
simple!

If you have multiple contiguous text files you want to put in the same hdf5 file, you can do it!
Simply pass all text files in the order you want them put into the hdf5 file. My code won't check
that they are contiguous (or consecutive), so be careful! 

Example usage:
    python convert_lammps_to_hdf5.py -i pos_1.dat pos_2.dat -o pos.hdf5
"""

import h5py, os, sys
import numpy as np
from itertools import islice
from timeit import default_timer
import argparse

# --------------------------------------------------------------------------------------------------

class c_timer:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,label,units='s'):

        """
        small tool for timing and printing timing info (copied from some other code of mine...)
        """ 

        self.label = label
        if units == 'm':
            self.units = 'm'
            self.scale = 1/60
        elif units == 'ms':
            self.units = 'ms'
            self.scale = 1000
        else:
            self.units = 's'
            self.scale = 1
        self.start_time = default_timer()

    # ----------------------------------------------------------------------------------------------

    def stop(self):

        """
        stop timer and print timing info
        """

        elapsed_time = default_timer()-self.start_time
        elapsed_time *= self.scale
        msg = f'\ntiming:   {self.label} {elapsed_time:9.5f} [{self.units}]'
        print(msg)

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

class c_compressor:
    
    # ----------------------------------------------------------------------------------------------
    
    def __init__(self,lammps_traj_files,hdf5_traj_file):
        
        """
        lammps_traj_files can be a list of files or a single one. all files are read and appended 
        to the same hdf5 file. if multiple are given, they should be consecutive and contain 
        exactly the same number of atoms in the same order. it doesnt matter how many steps are
        in each file, the code will check and work accordingly. 

        hdf5_traj_file is the output hdf5 file. if it exists, it is overwritten.
        """

        # i.e. length of data block that isnt the positions
        self.header_length = 9

        # output hdf5 file
        self.hdf5_traj_file = hdf5_traj_file

        # enforce list of input files is list (even if single file)
        if not isinstance(lammps_traj_files,list):
            lammps_traj_files = [lammps_traj_files]
        self.lammps_traj_files = lammps_traj_files
        self.num_lammps_files = len(self.lammps_traj_files)
        
        # check that the lammps files exist
        for _f in self.lammps_traj_files:
            if not os.path.exists(_f):
                print(f'\nthe file {_f} is missing!')
                sys.exit()

        print('\nlammps traj files:')
        print(self.lammps_traj_files)

        # get number of steps in each file... takes a while to 'loop' over the whole file
        self._get_meta_data()

        # now go and read the files/write to hdf5
        self._parse_text_files()

    # ----------------------------------------------------------------------------------------------

    def _convert_box_vecs(self,lat_params):

        """
        converts from lammps restricted triclinic box to generic lattice vectors.
        see https://docs.lammps.org/Howto_triclinic.html 

        xlo_bound xhi_bound xy
        ylo_bound yhi_bound xz
        zlo_bound zhi_bound yz

        take lat_params = [xlo xhi xy]
                          [ylo yhi xz]
                          [zlo zhi yz]
        and convert to 
               box_vecs = [ax   0   0]
                          [bx  by   0]
                          [cx  cy  cz]
        """

        _a = lat_params[0,1]-lat_params[0,0]
        _b = lat_params[1,1]-lat_params[1,0]
        _c = lat_params[2,1]-lat_params[2,0]

        _xy = lat_params[0,2]; _xz = lat_params[1,2]; _yz = lat_params[2,2]
        _cos_a = (_xy*_xz+_b*_yz)/(_b*_c); _cos_b = _xz/_c; _cos_g = _xy/_b

        _ax = _a
        _bx = _b*_cos_g
        _by = np.sqrt(_b**2-_bx**2)
        _cx = _c*_cos_b
        _cy = (_b*_c*_cos_a-_bx*_cx)/_by
        _cz = np.sqrt(_c**2-_cx**2-_cy**2)

        box_vecs = np.zeros((3,3),dtype=float)
        box_vecs[0,0] = _a
        box_vecs[1,0] = _bx; box_vecs[1,1] = _by
        box_vecs[2,0] = _cx; box_vecs[2,1] = _cy; box_vecs[2,2] = _cz

        return box_vecs

    # ----------------------------------------------------------------------------------------------

    def _process_types(self,types):

        """
        my code expects the types to start at 1 and be contiguous. this converts them to 
        contiguous and starting at 1
        """

        _unique, inverse = np.unique(types,return_inverse=True)
        return inverse+1

    # ----------------------------------------------------------------------------------------------

    def _parse_text_files(self):

        """
        loop over the text files, read each step one after another, and write them to the
        hdf5 database
        """

        _timer = c_timer('parse_text_files')

        # number of steps in each block to read (header+pos)
        _num_read = self.header_length+self.num_atoms

        # have types been written yet (only done at 1st step)
        _write_types = True

        # id, type, xu, yu, zu
        _data = np.zeros((self.num_atoms,5),dtype=float)
        _lat_params = np.zeros((3,3),dtype=float)

        with h5py.File(self.hdf5_traj_file,'w') as db:

            # initialize datasets
            _db_pos = db.create_dataset('cartesian_pos',
                                        shape=(self.num_steps,self.num_atoms,3),dtype=float)
            _db_types = db.create_dataset('types',shape=(self.num_atoms,),dtype=int)
            _db_box_vecs = db.create_dataset('box_vectors',shape=(self.num_steps,3,3),dtype=float)
            
            _step = 0

            # loop over the lammps files
            for ii, _lammps_file in enumerate(self.lammps_traj_files):
                print('\nnow reading file:',_lammps_file)
                print('this may take a while ...')

                # number of steps in this file
                _num_steps = self.num_steps_per_file[ii]

                # timer to read this file
                _f_timer = c_timer(f'read[{ii}]')
                
                with open(_lammps_file,'r') as _f:
                    for jj in range(_num_steps):

                        if jj % 100 == 0:
                            print(f'now on step {jj}/{_num_steps}')

                        # read this step
                        _lines = list(islice(_f,_num_read))

                        # lattice parameters
                        _lat_params[0,:] = _lines[5].strip().split()
                        _lat_params[1,:] = _lines[6].strip().split()
                        _lat_params[2,:] = _lines[7].strip().split()
                        _box_vecs = self._convert_box_vecs(_lat_params)

                        # convert types, positions, etc. to numbers
                        _lines = [_.strip().split() for _ in _lines[9:]]
                        _data[...] = _lines

                        # write to hdf5 file
                        _db_pos[_step,...] = _data[:,2:]
                        _db_box_vecs[_step,...] = _box_vecs[...]

                        # only write at the 1st step
                        if _write_types: 
                            _types = self._process_types(_data[:,1])
                            _db_types[...] = _types
                            _write_types = False

                        _step += 1

                _f_timer.stop()
            _timer.stop()
                        
    # ----------------------------------------------------------------------------------------------

    def _get_meta_data(self):

        """
        get number of atoms and number of steps from the lammps files
        """

        _timer = c_timer('get_meta_data')

        print('\ngetting meta_data from files ...')

        # make sure all files have the same number of atoms
        _num_atom_arr = np.zeros(self.num_lammps_files,dtype=int)

        self.num_lines_per_file = np.zeros(self.num_lammps_files,dtype=int)
        self.num_steps_per_file = np.zeros(self.num_lammps_files,dtype=int)
        
        for ii, _lammps_file in enumerate(self.lammps_traj_files):

            # get number of atoms in the file
            with open(_lammps_file,'r') as _f:
                for jj in range(3):
                    _f.readline()
                _num_atoms = int(_f.readline().strip())
            _num_atom_arr[ii] = _num_atoms
            
            # now check total number of lines/number of steps in the file
            with open(_lammps_file, "rb") as _f:
                _num_lines = sum(1 for _ in _f)
                _num_steps = int(_num_lines/(self.header_length+_num_atoms))

            self.num_steps_per_file[ii] = _num_steps
            self.num_lines_per_file[ii] = _num_lines
        
            print('\nfile name:',_lammps_file)
            print('num atoms:',_num_atoms)
            print('num lines:',_num_lines)
            print('num steps:',_num_steps)

        # make sure all files have same number of atoms
        self.num_atoms = _num_atom_arr[0]
        if any(self.num_atoms-_num_atom_arr):
            print('\nnumber of atoms in all files must be the same!')
            sys.exit()

        self.num_steps = self.num_steps_per_file.sum()
        print('\ntotal num steps:',self.num_steps)

        _timer.stop()

    # ----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

def parse_args():

    """
    get lammps and hdf5 file names from command line
    """

    # get cmd line args
    description = 'command line args for \'compress.py\''
    cmd_parser = argparse.ArgumentParser(description=description)

    # input file
    help_msg = 'lammps trajectory files to be merged into an hdf5 database. give one or\n'     \
               'several. they should be consecutive (in time) and should have exactly\n'       \
               'the same number of atoms in the same order. it doenst matter how many steps\n' \
               'are in each however\nthey should have be written with a command like:\n'       \
               '---- dump          pos all custom ${dump_freq} pos.dat id type xu yu zu\n'      \
               '---- dump_modify   pos sort id'
    cmd_parser.add_argument('-i','--input-files',default=['pos.dat'],help=help_msg,nargs='+')
    
    # output file
    help_msg = 'the name of the output hdf5 file'
    cmd_parser.add_argument('-o','--output-file',default='pos.hdf5',help=help_msg)

    # get cmd line args
    cmd_args = cmd_parser.parse_args()
    input_files = cmd_args.input_files
    output_file = cmd_args.output_file

    return input_files, output_file

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    lammps_traj_files, hdf5_traj_file = parse_args()
    c_compressor(lammps_traj_files,hdf5_traj_file)













