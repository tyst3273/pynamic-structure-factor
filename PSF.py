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

"""

# system modules
import argparse

# custom modules
from psf.m_PSF import c_PSF

# --------------------------------------------------------------------------------------------------

# get cmd line args
description = 'command line args for \'PSF.py\''
cmd_parser = argparse.ArgumentParser(description=description)

# input file
help_msg = 'input file for \'PSF.py\'. it is a python file that is imported as a module'
cmd_parser.add_argument('-i','--input-file',default='psf_input.py',help=help_msg)

# get cmd line args
cmd_args = cmd_parser.parse_args()
input_file = cmd_args.input_file

# --------------------------------------------------------------------------------------------------

PSF = c_PSF(input_file)
PSF.standard_run()

# --------------------------------------------------------------------------------------------------


