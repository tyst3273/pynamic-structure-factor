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
perform a 'standard' calculation by reading config file (file path is read from 
cmd-line or defaults inside of c_config), doing what is requested, then writing 
the results and exiting.

run it as e.g.:
    python PSF.py

or:
    python PSF.py -i input_params.py

or: 
    python PSF.py -i input_1.py input_2.py ... 

"""

# custom modules
from psf.m_PSF import c_PSF
from psf.m_config import get_input_files

# --------------------------------------------------------------------------------------------------

# get file(s) from cmd line
input_files = get_input_files()

# loop over files 
for input_file in input_files:

    # set up PSF API
    PSF = c_PSF(input_file)

    # do standard run from file
    PSF.standard_run()

# --------------------------------------------------------------------------------------------------


