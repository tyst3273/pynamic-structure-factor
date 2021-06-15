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

mod_path = 'modules'
sys.path.append(mod_path)

import mod_invars
import mod_io
import mod_lattice
import mod_Qpoints
import mod_sqw
import mod_utils 

# -----------------------------------------------------------------------

# get the input file from command line
if len(sys.argv) != 1:
    input_files = sys.argv[1:]
else:
    input_files = ['input_params']

# -----------------------------------------------------------------------

# stub for mpi
rank = 0
num_ranks = 1

# -----------------------------------------------------------------------

# print startup info
if rank == 0:

    # start the outer timer
    outer_start = timer()

    # print herald
    mod_utils.print_herald(num_ranks)

    # print the number of files
    num_files = len(input_files)
    message = f'there are {num_files} files to do:\n'
    for input_file in input_files:
        message = message+f' - {input_file}\n'  
    mod_io.print_stdout(message,msg_type='FILES')

# -----------------------------------------------------------------------

# loop over input files
for input_file in input_files:

    # only rank 0 needs to read file, run timer, etc
    if rank == 0:

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

        # setup Qpoints or read the file
        Qpoints = mod_Qpoints.Qpoints()
        Qpoints.generate_Qpoints(invars)

        # split Qpoints over procs
        Qpoints.distribute_Q_over_procs(invars,num_ranks)

        # distribute the vars over all procs
        for ii in range(1,num_ranks):
            comm.send(invars,dest=ii,tag=0)
            comm.send(lattice,dest=ii,tag=1)
            comm.send(Qpoints,dest=ii,tag=2)
    
    # ----------------------------------------------------------------

    # other ranks receive data produced on rank 0
    else:

        # receive the vars from rank 0
        invars = comm.recv(source=0,tag=0)
        lattice = comm.recv(source=0,tag=1)
        Qpoints = comm.recv(source=0,tag=2)

    # ----------------------------------------------------------------

    # get Qpoints that are on this rank, convert to 1/A
    Qpoints.rank_init(lattice,rank)

    # open the hdf5 file
    traj_file = mod_io.traj_file(invars)

    # setup calculator
    sqw = mod_sqw.sqw(invars,Qpoints,rank)

    # run the calculation
    sqw.calculate(invars,Qpoints,lattice,traj_file)

    # close the hdf5 files
    traj_file.close()

    # ----------------------------------------------------------------

    # rank 0 gathers all the results to 1 array and writes it
    if rank == 0:

        # gather all of the SQW results
        if invars.compute_sqw:

            # get from all ranks
            sqw_from_ranks = [sqw.sqw]
            for ii in range(1,num_ranks):
                sqw_ii = comm.recv(source=ii,tag=11)
                sqw_from_ranks.append(sqw_ii)

            # assemble the SQW results back into one
            sqw_total = np.zeros((sqw.num_freq,Qpoints.total_Qsteps))
            mod_utils.assemble_sqw(sqw_from_ranks,sqw_total)

            # save total SQW to final hdf5 file
            f_name = invars.outfile_prefix+f'_SQW_FINAL.hdf5'
            mod_io.save_sqw(invars,Qpoints.total_reduced_Q,sqw.meV,sqw_total,f_name)

        # ---------------------------------------------------------------------

        # gather and save bragg results
        if invars.compute_bragg:
            
            # get from all ranks
            bragg_from_ranks = [sqw.bragg]
            for ii in range(1,num_ranks):
                bragg_ii = comm.recv(source=ii,tag=12)
                bragg_from_ranks.append(bragg_ii)

            # assemble into a single array
            bragg_total = np.zeros(Qpoints.total_Qsteps)
            mod_utils.assemble_timeavg(bragg_from_ranks,bragg_total)

            # write to final file
            f_name = invars.outfile_prefix+f'_BRAGG_FINAL.hdf5'
            mod_io.save_bragg(invars,Qpoints.total_reduced_Q,bragg_total,f_name)

        # ---------------------------------------------------------------------

        # gather and save timeavg results
        if invars.compute_timeavg:

            # get from all ranks
            timeavg_from_ranks = [sqw.timeavg]
            for ii in range(1,num_ranks):
                timeavg_ii = comm.recv(source=ii,tag=13)
                timeavg_from_ranks.append(timeavg_ii)

            # assemble into a single array
            timeavg_total = np.zeros(Qpoints.total_Qsteps)
            mod_utils.assemble_timeavg(timeavg_from_ranks,timeavg_total)

            # write to final file
            f_name = invars.outfile_prefix+f'_TIMEAVG_FINAL.hdf5'
            mod_io.save_timeavg(invars,Qpoints.total_reduced_Q,timeavg_total,f_name)

        # -------------------------------------------------------------------------

        # calculate and print elapsed time for file
        inner_end = timer()
        inner_time = (inner_end-inner_start)/60 #minutes
        message = f'elapsed time for this file: {inner_time:2.3f} minutes'
        mod_io.print_stdout(message,msg_type='TIMING')

    # ---------------------------------------------------------------

    # other ranks send thier data to rank 0
    else:

        # send the SQW results to rank 0
        if invars.compute_sqw:
            comm.send(sqw.sqw,dest=0,tag=11)

        # send the bragg intensity too
        if invars.compute_bragg:
            comm.send(sqw.bragg,dest=0,tag=12)

        # send the timeavg intensity too
        if invars.compute_timeavg:
            comm.send(sqw.timeavg,dest=0,tag=13)

    # ----------------------------------------------------------------

# print total elapsed time
if rank == 0:

    # calculate and print total elapsed time
    outer_end = timer()
    outer_time = (outer_end-outer_start)/60 #minutes
    message = f'total elapsed time: {outer_time:2.3f} minutes'
    mod_io.print_stdout(message,msg_type='TIMING')

# -----------------------------------------------------------------------------------









