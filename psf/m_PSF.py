#!/home/ty/anaconda3/bin/python3

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

import psf.m_communicator as m_communicator
import psf.m_config as m_config
import psf.m_timing as m_timing
import psf.m_io as m_io
from psf.m_printing import print_preamble, print_goodbye


class c_PSF:

    # ----------------------------------------------------------------------------------------------
    
    def __init__(self,input_file=None):

        """
        main class that holds 'macros' to do stuff
        """
        
        # print author info etc.
        print_preamble()

        # timers
        self.timers = m_timing.c_timers()

        # options for the calculation
        self.config = m_config.c_config(input_file)

    # ----------------------------------------------------------------------------------------------

    def setup_calculation(self,pos=None,types=None,**kwargs):

        """
        set up the 'communicator' object to pass stuff back and forth
        """

        # read input file and overwrite input args in file with kwargs
        self.config.set_config(**kwargs)

        # 'communicator' to conveniently pass objects and data around
        self.comm = m_communicator.c_communicator(self.config,self.timers)
        self.comm.setup_calculation(pos,types)

    # ----------------------------------------------------------------------------------------------

    def standard_run(self):

        """
        the usual way to run the code. it reads the config file, does what was requested, 
        writes the results, the exits.
        """

        # read input file and setup communicator
        self.setup_calculation()

        # run the calculation 
        self.run()

    # ----------------------------------------------------------------------------------------------

    def finalize(self):

        """
        set up the 'communicator' object to pass stuff back and forth
        """

        # write output files
        self.write_strufacs()
            
        # print timing
        self.timers.print_timing()

        # print 'goodbye' message
        print_goodbye()

    # ----------------------------------------------------------------------------------------------

    def write_strufacs(self):

        """
        write the data to files
        """

        # setup io object to write files
        writer = m_io.c_writer(self.config,self.comm)
        writer.write_structure_factors()
    
    # ----------------------------------------------------------------------------------------------

    def run(self):

        """
        run the calculation
        """

        self.timers.start_timer('PSF',units='m')

        # this loops over all blocks, calculating errything
        self.comm.strufacs.calculate_structure_factors()

        # if calculated on a mesh, unfold mesh onto full reciprocal space
        if self.comm.Qpoints.use_mesh:
            self.comm.strufacs.put_on_mesh()

        self.timers.stop_timer('PSF')

        self.finalize()

    # ----------------------------------------------------------------------------------------------





