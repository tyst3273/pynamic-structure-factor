import numpy as np
from FileIO import print_stdout


# =======================================================================
# ------------ subroutines to split Q points over processes. -------------
# =======================================================================

def prepare_Qpoints(params,n_ranks):

    """
    num_Q_per_proc is determined as the largest integer diving the total_Qsteps number. the remainder
    is placed on rank 0 (if there is a remainder...)
    """

    total_reduced_Q = params.total_reduced_Q 
    total_Qsteps = params.total_Qsteps
    
    # if only 1 proc, return the total Q set
    if n_ranks == 1:
        message = f'process: 0    Q points: {total_Qsteps:g}\n'
        print_stdout(message,msg_type='Q points on each process')
        return [total_reduced_Q]

    # else, take num_Q_per_proc to be the largest integer dividing total_Qsteps
    num_Q_per_proc = total_Qsteps//n_ranks

    nQpp_arr = np.zeros(n_ranks).astype(int)
    nQpp_remainder = total_Qsteps-n_ranks*num_Q_per_proc 

    nQpp_arr[:] = num_Q_per_proc
    nQpp_arr[0] = nQpp_arr[0]+nQpp_remainder

    if nQpp_remainder != 0:
        message = (f'one proc. will compute {nQpp_remainder+num_Q_per_proc} Qpoints;'
                    f' the rest will compute {num_Q_per_proc}.\n this might be inefficient.')
        print_stdout(message,msg_type='WARNING (PARALLELISM)')

    else:
        message = f'each proc. will compute {num_Q_per_proc} Qpoints'
        print_stdout(message,msg_type='NOTE (PARALLELISM)')

    message = f'process: 0    Q points: {nQpp_arr[0]:g}\n'
    for ii in range(1,n_ranks):
        message = message+f' process: {ii:g}    Q points: {nQpp_arr[ii]:g}\n'
    print_stdout(message,msg_type='Q points on each process')

    shift = 0
    reduced_Q_set = []
    for ii in range(n_ranks):
        reduced_Q_set.append(total_reduced_Q[shift:shift+nQpp_arr[ii],:])
        shift = shift+nQpp_arr[ii]

    return reduced_Q_set


# =======================================================================
# ------- subroutine to combine SQW from each proc into one -------------
# =======================================================================

def assemble_sqw_total(params,sqw_total_set,n_ranks):

    """
    assemble the sqw data from each proc into a single array
    """

    sqw_total = np.zeros((params.num_freq,params.total_Qsteps))

    shift = 0
    for ii in range(n_ranks):
        nQpp = sqw_total_set[ii].shape[1]
        sqw_total[:,shift:shift+nQpp] = sqw_total_set[ii]
        shift = shift+nQpp

    return sqw_total



# =======================================================================
# ------------ cleanest way i could figure out how to exit --------------
# =======================================================================

class ParalExcept(Exception):

    """
    to exit the code, just run 'raise ParalExcept'. note, must call python with
    python -m mpi4py [script] or it will only kill the calling process, possibly
    resulting in deadlock.
    """

    def __init__(self,message='an error occured. aborting',rank=0):

        if rank == 0:
            print_stdout(message,msg_type='ERROR')











