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

"""
'driver' script for psf code. gonna put 'tasks' in here, i.e. macros for calculating 
all kinds of stuff in here... but for now, its just the main script.
"""

import psf.m_communicator as m_communicator
import psf.m_config as m_config
import psf.m_timing as m_timing


class c_PSF:

    # ----------------------------------------------------------------------------------------------
    
    def __init__(self,input_file=None):

        """
        main class that holds 'macros' to do stuff
        """

        self.print_preamble()

        # timers
        self.timers = m_timing.c_timers()

        # options for the calculation
        self.config = m_config.c_config(input_file)

    # ----------------------------------------------------------------------------------------------

    def print_preamble(self):

        """
        self explanatory ...
        """

        preamble = '\n\n#######################################################################\n'
        preamble += 'Pynamic Structure Factors \n'
        preamble += 'author: Tyler C. Sterling\n'
        preamble += 'email: ty.sterling@colorado.edu\n'
        preamble += 'affil: Physics Dept., University of Colorado Boulder\n'
        preamble += '  Neurtron Scattering and Raman Spectroscopy Lab\n'
        preamble += '#######################################################################\n'
        print(preamble)

    # ----------------------------------------------------------------------------------------------

    def print_goodbye(self):

        """
        self explanatory ...
        """

        self.timers.print_timing()

        goodbye = '#######################################################################\n'
        goodbye += 'the calculation finished willy-nilly\n' 
        goodbye += 'as always, check the results carefully\n'
        goodbye += 'bye!\n'
        goodbye += '#######################################################################\n'
        print(goodbye)

    # ----------------------------------------------------------------------------------------------

    def setup_communicator(self,pos=None,types=None):

        """
        set up the 'communicator' object to pass stuff back and forth
        """

        # 'communicator' to conveniently pass objects in/out of stuff
        self.comm = m_communicator.c_communicator(self.config,self.timers)

        if self.config.trajectory_format in ['external']:
            self.comm.setup_calculation(pos,types)
        else:
            self.comm.setup_calculation()

    # ----------------------------------------------------------------------------------------------

    def standard_run(self):

        """
        the usual way to run the code. it reads the config file, does what was requested, 
        writes the results, the exits.
        """

        self.timers.start_timer('standard_run',units='m')

        # read input file 
        self.get_config_from_file()

        # 'communicator' to conveniently pass objects in/out of stuff
        self.setup_communicator()

        # run the calculation 
        self.run()

        # write output files
        self.write_strufacs()

        self.timers.stop_timer('standard_run')

    # ----------------------------------------------------------------------------------------------

    def write_strufacs(self):

        """
        write the data to files
        """

        # write output files
        self.comm.writer.write_structure_factors()
    
    # ----------------------------------------------------------------------------------------------

    def run(self):

        """
        run the calculation
        """

        # this loops over all blocks, calculating errything
        self.comm.strufacs.calculate_structure_factors()

        # unfold mesh onto full reciprocal space
        if self.comm.Qpoints.use_mesh:
            self.comm.strufacs.put_on_mesh()

    # ----------------------------------------------------------------------------------------------

    def get_config_from_file(self):

        """
        get config info from file
        """

        self.config.set_config()

    # ----------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    """
    if not being imported as a module, perform a 'standard' calculation by reading config file
    (file path is read from cmd-line or defaults in c_config), doing what is requested, then
    writing the results and exiting
    """

    PSF = c_PSF()
    PSF.standard_run()
    PSF.print_goodbye()

# --------------------------------------------------------------------------------------------------




