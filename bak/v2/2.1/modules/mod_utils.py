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

# -----------------------------------------------------------------------------------------

class PSF_exception(Exception):
    """
    to exit the code, just run 'raise ParalExcept'. note, must call python with
    python -m mpi4py [script] or it will only kill the calling process, possibly
    resulting in deadlock if during a send/recv
    """
    def __init__(self,message='an error occured. aborting',rank=0):
        if rank == 0:
            print_stdout(message,msg_type='ERROR')

# -----------------------------------------------------------------------------------------

def print_stdout(message,msg_type=None):
    """
    print to screen using common formatting
    """
    if msg_type != None:
        print(f'\n ** {msg_type} **')
    print(f' {message}',flush=True)

# -----------------------------------------------------------------------------------------

def print_herald(num_ranks):
    """
    print author info etc to screen at startup
    """
    herald = """
 Pynamic Structure Factor, version 2.1

 Now with MPI support!

 Author: Ty Sterling
         Department of Physics
         University of Colorado Boulder

 Email: ty.sterling@colorado.edu

 """

    banner = """
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""
    print(herald)
    print(f' ** running the S(Q,w) calculation using {num_ranks} processes **\n')
    print(banner,flush=True)


# -----------------------------------------------------------------------------------------

def assemble_sqw(sqw_from_ranks,sqw_total):
    """
    assemble SQW from all the procs back into 1 array.
    """
    shift = 0
    num_ranks = len(sqw_from_ranks)
    for ii in range(num_ranks):
        nQpp = sqw_from_ranks[ii].shape[1]
        sqw_total[:,shift:shift+nQpp] = sqw_from_ranks[ii]
        shift = shift+nQpp

# -----------------------------------------------------------------------------------------

def assemble_timeavg(timeavg_from_ranks,timeavg_total):
    """
    assemble timeaveraged (or bragg) intensity from all the procs back into 1 array.
    """
    shift = 0
    num_ranks = len(timeavg_from_ranks)
    for ii in range(num_ranks):
        nQpp = timeavg_from_ranks[ii].shape[0]
        timeavg_total[shift:shift+nQpp] = timeavg_from_ranks[ii]
        shift = shift+nQpp

# ------------------------------------------------------------------------------------------






