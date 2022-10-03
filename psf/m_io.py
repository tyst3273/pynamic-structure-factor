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


class c_writer:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm):

        """
        write file containing structure factors
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config

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

        msg = f'\n*** io ***\nwriting data to\n  \'{self.output_file}\''
        print(msg)
        
        # write header now; otherwise the other methods do if called independently
        self._write_header()

        if self.config.calc_sqw:
            self.write_sqw()
        if self.config.calc_diffuse:
            self.write_diffuse()
        if self.config.calc_bragg:
            self.write_bragg()
        if self.config.calc_rho_squared:
            self.write_rho_squared()

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

                db.create_dataset('lattice_vectors',shape=(3,3),dtype=float)
                db['lattice_vectors'][...] = self.comm.lattice.lattice_vectors[...]

                db.create_dataset('reciprocal_lattice_vectors',shape=(3,3),dtype=float)
                db['reciprocal_lattice_vectors'][...] = \
                        self.comm.lattice.reciprocal_lattice_vectors[...]
    
                db.create_dataset('mesh',shape=(1),dtype=int)

                if self.comm.Qpoints.use_mesh:
                    db['mesh'][...] = int(True)
                    db.create_dataset('H',shape=_Qpoints.num_H,dtype=float)
                    db['H'][...] = _Qpoints.H[...]
                    db.create_dataset('K',shape=_Qpoints.num_K,dtype=float)
                    db['K'][...] = _Qpoints.K[...]
                    db.create_dataset('L',shape=_Qpoints.num_L,dtype=float)
                    db['L'][...] = _Qpoints.L[...]
                else:
                    db['mesh'][...] = int(False)
                    db.create_dataset('Q_rlu',shape=_Qpoints.Q_rlu.shape,dtype=float)
                    db['Q_rlu'][...] = _Qpoints.Q_rlu[...]
                    db.create_dataset('Q_cart',shape=_Qpoints.Q_cart.shape,dtype=float)
                    db['Q_cart'][...] = _Qpoints.Q_cart[...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------
    
    def write_sqw(self):

        """
        write S(Q,w) info to file
        """

        self._write_header()

        try:
            with h5py.File(self.output_file,'a') as db:

                db.create_dataset('sqw',shape=self.comm.strufacs.sqw.shape,dtype=float)
                db['sqw'][...] = self.comm.strufacs.sqw[...]

                db.create_dataset('energy',shape=self.comm.strufacs.energy.shape,dtype=float)
                db['energy'][...] = self.comm.strufacs.energy[...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def write_bragg(self):

        """
        write bragg info to file
        """

        self._write_header()

        try:
            with h5py.File(self.output_file,'a') as db:

                db.create_dataset('sq_bragg',shape=self.comm.strufacs.sq_bragg.shape,dtype=float)
                db['sq_bragg'][...] = self.comm.strufacs.sq_bragg[...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def write_rho_squared(self):

        """
        write |exp(iQ.r(t))|**2 to file
        """

        self._write_header()

        try:
            with h5py.File(self.output_file,'a') as db:

                db.create_dataset('rho_sq', \
                        shape=self.comm.strufacs.rho_sq.shape,dtype=float)
                db['rho_sq'][...] = self.comm.strufacs.rho_sq[...]

                db.create_dataset('time',self.comm.traj.time.shape,dtype=float)
                db['time'][...] = self.comm.traj.time[...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def write_diffuse(self):

        """
        write S(Q,w) info to file
        """

        self._write_header()

        try:
            with h5py.File(self.output_file,'a') as db:

                db.create_dataset('sq_diffuse', \
                        shape=self.comm.strufacs.sq_diffuse.shape,dtype=float)
                db['sq_diffuse'][...] = self.comm.strufacs.sq_diffuse[...]

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

    def read_bragg(self):

        """
        read bragg intensity
        """

        self._read_header()

        try:
            with h5py.File(self.input_file,'r') as db:
                self.sq_bragg = db['sq_bragg'][...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def read_diffuse(self):

        """
        read diffuse intensity
        """

        self._read_header()

        try:
            with h5py.File(self.input_file,'r') as db:
                self.sq_diffuse = db['sq_diffuse'][...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------

    def read_rho_squared(self):

        """
        read rho squared
        """

        self._read_header()

        try:
            with h5py.File(self.input_file,'r') as db:
                self.rho_sq = db['rho_sq'][...]
                self.time = db['time'][...]

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------






