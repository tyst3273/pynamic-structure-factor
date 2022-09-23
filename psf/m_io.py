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

from psf.m_error import crash


class c_io:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm):

        """
        read/write file containing structure factors
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config

        # whether or no 'header' has been written
        self.header_flag = False

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

        msg = f'\n*** io ***\nwriting data to\n  \'{self.output_file}\''
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def write_structure_factors(self):

        """
        call methods to write everything
        """

        if self.config.calc_sqw:
            self.write_sqw()
        if self.config.calc_diffuse:
            self.write_diffuse()
        if self.config.calc_bragg:
            self.write_bragg()

    # ----------------------------------------------------------------------------------------------

    def _write_header(self):

        """
        write S(Q,w) info to file
        """

        if self.header_flag:
            return

        self.header_flag = True

        try:
            with h5py.File(self.output_file,'a') as db:

                if self.comm.Qpoints.use_mesh:
                    db['mesh'] = True
                else:
                    db['mesh'] = False

        except Exception as _ex:
            crash(self.crash_msg,_ex)

    # ----------------------------------------------------------------------------------------------
    
    def write_sqw(self):

        """
        write S(Q,w) info to file
        """

        self._write_header()

    # ----------------------------------------------------------------------------------------------

    def read_sqw(self):

        """
        read S(Q,w) info from file
        """

        pass

    # ----------------------------------------------------------------------------------------------

    def write_bragg(self):

        """
        write bragg info to file
        """

        self._write_header()

    # ----------------------------------------------------------------------------------------------

    def read_bragg(self):

        """
        read bragg info from file
        """

        pass

    # ----------------------------------------------------------------------------------------------

    def write_diffuse(self):

        """
        write S(Q,w) info to file
        """

        self._write_header()

    # ----------------------------------------------------------------------------------------------

    def read_diffuse(self):

        """
        read S(Q,w) info from file
        """

        pass

    # ----------------------------------------------------------------------------------------------

