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

from timeit import default_timer as timer
import numpy as np
import sys

import modules.mod_invars as mod_invars
import modules.mod_io as mod_io
import modules.mod_lattice as mod_lattice
import modules.mod_Qpoints as mod_Qpoints
import modules.mod_sqw as mod_sqw
import modules.mod_utils  as mod_utils

# --------------------------------------------------------------------------------

# get the input file from command line
if len(sys.argv) != 1:
    input_files = sys.argv[1:]
else:
    input_files = ['input_params']

# --------------------------------------------------------------------------------

# start the outer timer
outer_start = timer()

# print herald
mod_utils.print_herald()

# print the number of files
num_files = len(input_files)
message = f'there are {num_files} files to do:\n'
for input_file in input_files:
    message = message+f' - {input_file}\n'  
mod_io.print_stdout(message,msg_type='FILES')

# --------------------------------------------------------------------------------

# loop over input files
for input_file in input_files:

    # print which file
    message = '\n------------------------------------------'
    mod_io.print_stdout(message)
    message = f' now working on the file \'{input_file}\''
    mod_io.print_stdout(message,msg_type='FILES')

    # start a timer
    inner_start = timer()

    # read the input file
    invars = mod_invars.input_variables()
    invars.parse_input(input_file)

    # setup the lattice
    lattice = mod_lattice.lattice(invars)

    # setup Qpoints along path or read the file
    Qpoints = mod_Qpoints.Qpoints()
    Qpoints.generate_Qpoints(invars,lattice)

    # distribute Qpoints over procs for multiprocessing
    Qpoints.distribute_Q_over_procs(invars.num_processes)

    # open the hdf5 file
    traj_file = mod_io.traj_file(invars)

    # setup calculator
    sqw = mod_sqw.sqw(invars,Qpoints)

    # ** run the calculation **
    sqw.calculate(invars,Qpoints,lattice,traj_file)

    # close the hdf5 files
    traj_file.close()

    # ----------------------------------------------------------------------------

    # save the SQW results
    if invars.compute_sqw:

        # save total SQW to final hdf5 file
        f_name = invars.outfile_prefix+f'_SQW_FINAL.hdf5'
        mod_io.save_sqw(invars,Qpoints.total_reduced_Q,sqw.meV,sqw.sqw,f_name)

    # ----------------------------------------------------------------------------

    # save bragg results
    if invars.compute_bragg:
            
        # write to final file
        f_name = invars.outfile_prefix+f'_BRAGG_FINAL.hdf5'
        mod_io.save_bragg(invars,Qpoints.total_reduced_Q,sqw.bragg,f_name)

    # ----------------------------------------------------------------------------

    # save timeavg results
    if invars.compute_timeavg:

        # write to final file
        f_name = invars.outfile_prefix+f'_TIMEAVG_FINAL.hdf5'
        mod_io.save_timeavg(invars,Qpoints.total_reduced_Q,sqw.timeavg,f_name)

    # ----------------------------------------------------------------------------

    # calculate and print elapsed time to do this file
    inner_end = timer()
    inner_time = (inner_end-inner_start)/60 # minutes
    message = f'elapsed time for this file: {inner_time:2.3f} minutes'
    mod_io.print_stdout(message,msg_type='TIMING')

# --------------------------------------------------------------------------------

# calculate and print total elapsed time
outer_end = timer()
outer_time = (outer_end-outer_start)/60 # minutes
message = f'total elapsed time: {outer_time:2.3f} minutes'
mod_io.print_stdout(message,msg_type='TIMING')

# --------------------------------------------------------------------------------









