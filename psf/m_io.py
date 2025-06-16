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

import os
import numpy as np
import h5py 

from psf.m_error import crash, check_file

# --------------------------------------------------------------------------------------------------

class c_writer:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm):

        """
        write file containing structure factors
        """
        
        # copy refs to stuff
        self.comm = comm
        self.config = config

        # wheter or not to write output file at all
        self.write_output_file = True
        if self.config.output_prefix is None:
            self.write_output_file = False
            return

        # whether or not 'header' has been written to file
        self.write_flag = False

        # get the file name
        self.output_prefix = self.config.output_prefix+'_'
        self.output_directory = self.config.output_directory
        self.output_file = self.output_prefix+'STRUFACS.hdf5'
        self.output_file = os.path.join(self.output_directory,self.output_file)

        self.crash_msg = 'could not write ouput file! check destination\n' \
                         'and h5py installation and try again!\n'

        # make directory if it doesnt exist
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)

        # remove file if it already exists
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

    # ----------------------------------------------------------------------------------------------

    def write_structure_factors(self):

        """
        call methods to write everything
        """

        if not self.write_output_file:
            msg = f'\n*** io ***\nnot writing output file'
            print(msg)
            return

        msg = f'\n*** io ***\nwriting data to\n  \'{self.output_file}\''
        print(msg)
        
        # write header now; otherwise the other methods do if called independently
        self._write_header()

        self._write_elastic()

        if self.config.calc_sqw:
            self._write_sqw()

    # ----------------------------------------------------------------------------------------------

    def _write_header(self):

        """
        write common info for all calcs
        """

        if self.write_flag:
            return

        self.write_flag = True

        _Qpoints = self.comm.Qpoints

        try:
            with h5py.File(self.output_file,'a') as db:

                db.create_dataset('lattice_vectors',
                        data=self.comm.lattice.lattice_vectors)

                db.create_dataset('reciprocal_lattice_vectors',
                        data=self.comm.lattice.reciprocal_lattice_vectors)
    
                db.create_dataset('mesh',shape=(1),dtype=int)

                if self.comm.Qpoints.use_mesh:
                    db['mesh'][...] = int(True)
                    db.create_dataset('H',data=_Qpoints.H)
                    db.create_dataset('K',data=_Qpoints.K)
                    db.create_dataset('L',data=_Qpoints.L)
                else:
                    db['mesh'][...] = int(False)
                    db.create_dataset('Q_rlu',data=_Qpoints.Q_rlu)
                    db.create_dataset('Q_cart',data=_Qpoints.Q_cart)

                if hasattr(_Qpoints,'Q_path_verts'):
                    db.create_dataset('Q_path_verts',data=_Qpoints.Q_path_verts)

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------
    
    def _write_sqw(self):

        """
        write S(Q,w) info to file
        """

        self._write_header()

        try:
            with h5py.File(self.output_file,'a') as db:
                db.create_dataset('sqw',data=self.comm.strufacs.sqw)
                db.create_dataset('energy',data=self.comm.strufacs.energy)

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def _write_elastic(self):

        """
        write elastic info to file
        """

        self._write_header()

        try:
            with h5py.File(self.output_file,'a') as db:
                db.create_dataset('sq_elastic',data=self.comm.strufacs.sq_elastic)

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------------


class c_reader:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,input_file='psf_STRUFACS.hdf5'):

        """
        read file containing structure factors
        """

        check_file(input_file)
        self.input_file = input_file

        self.read_flag = False

        self.crash_msg = 'could not write ouput file! check destination\n' \
                         'and h5py installation and try again!\n'

    # ----------------------------------------------------------------------------------------------

    def _read_header(self):

        """
        read common info from all calcs
        """

        if self.read_flag:
            return

        self.read_flag = True

        try:
            with h5py.File(self.input_file,'r') as db:

                self.lattice_vectors = db['lattice_vectors'][...]
                self.reciprocal_lattice_vectors = db['reciprocal_lattice_vectors'][...]

                if db['mesh'][...]:
                    self.mesh = True
                    self.H = db['H'][...]
                    self.K = db['K'][...]
                    self.L = db['L'][...]
                else:
                    self.mesh = False
                    self.Q_rlu = db['Q_rlu'][...]
                    self.Q_cart = db['Q_cart'][...]
                
                if 'Q_path_verts' in db.keys():
                    self.Q_path_verts = db['Q_path_verts'][...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def read_sqw(self):

        """
        read sqw 
        """

        self._read_header()

        try:
            with h5py.File(self.input_file,'r') as db:
                self.sqw = db['sqw'][...]
                self.energy = db['energy'][...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def read_elastic(self):

        """
        read elastic intensity
        """

        self._read_header()

        try:
            with h5py.File(self.input_file,'r') as db:
                self.sq_elastic = db['sq_elastic'][...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def get_energy_inds(self,e_min,e_max):
        
        """
        get the indices of sqw in specified energy range
        """

        inds = np.flatnonzero(self.energy <= e_max)
        inds = np.intersect1d(np.flatnonzero(self.energy >= e_min),inds)
        return inds

    # ----------------------------------------------------------------------------------------------

    def read_energy_integrated_sqw(self,e_min=None,e_max=None):

        """
        integrate (sum) sqw over specified energy range
        """

        self._read_header()

        try:
            with h5py.File(self.input_file,'r') as db:
                _sqw = db['sqw'][...]
                self.energy = db['energy'][...]
        except Exception as _ex:
            crash(self.crash_msg,_ex)

        if e_min is None:
            e_min = self.energy.min()-1
        if e_max is None:
            e_max = self.energy.max()+1

        e_inds = self.get_energy_inds(e_min,e_max)

        self.sq_integrated = _sqw[...,e_inds].sum(axis=-1)

    # ----------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------

def merge_strufac_files(strufac_files,output_file):

    """
    take a list of strufac files and merge them all together, averaging over sqw and sq_elastic 
    if present
    """

    if not isinstance(strufac_files,list):
        msg = 'strufac_files should be a list of file names'
        crash(msg)
        return

    for file in strufac_files:
        check_file(file)

    num_files = len(strufac_files)

    has_sqw = False
    init = True

    msg = '\nmerging structure factor files:'
    print(msg)

    # read all of the files and average sq_elastic and sqw
    for file in strufac_files:

        msg = f'  {file}'
        print(msg)    

        with h5py.File(file,'r') as db:

            # if reading the first file, initialize arrays
            if init:

                lattice_vectors = db['lattice_vectors'][...]
                reciprocal_lattice_vectors = db['reciprocal_lattice_vectors'][...]
                mesh = db['mesh'][...]

                sq_elastic = db['sq_elastic'][...]

                if 'sqw' in db.keys():
                    has_sqw = True
                    sqw = db['sqw'][...]
                    energy = db['energy'][...]

                if mesh == True:
                    H = db['H'][...]
                    K = db['K'][...]
                    L = db['L'][...]
                else:
                    Q_rlu = db['Q_rlu'][...]
                    Q_cart = db['Q_cart'][...]

                if 'Q_path_verts' in db.keys():
                    Q_path_verts = db['Q_path_verts'][...]

                init = False

            else:

                sq_elastic += db['sq_elastic'][...]

                if has_sqw:
                    sqw += db['sqw'][...]

    msg = f'\ninto file:\n  {output_file}\n'
    print(msg)

    # write the merged data to a new file  
    with h5py.File(output_file,'w') as db:

        db.create_dataset('lattice_vectors',data=lattice_vectors)
        db.create_dataset('reciprocal_lattice_vectors',data=reciprocal_lattice_vectors)

        if mesh:
            db.create_dataset('mesh',data=True)
            db.create_dataset('H',data=H)
            db.create_dataset('K',data=K)
            db.create_dataset('L',data=L)
        else:
            db.create_dataset('mesh',data=False)
            db.create_dataset('Q_rlu',data=Q_rlu)
            db.create_dataset('Q_cart',data=Q_cart)

        if 'Q_path_verts' in locals():
            db.create_dataset('Q_path_verts',data=Q_path_verts)

        sq_elastic /= num_files
        db.create_dataset('sq_elastic',data=sq_elastic)

        if has_sqw:
            sqw /= num_files
            db.create_dataset('sqw',data=sqw)
            db.create_dataset('energy',data=energy)

# --------------------------------------------------------------------------------------------------


