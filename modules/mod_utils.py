
# -----------------------------------------------------------------------------------------

class PSF_exception(Exception):
    """
    to exit the code, just run 'raise ParalExcept'. note, must call python with
    python -m mpi4py [script] or it will only kill the calling process, possibly
    resulting in deadlock.
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
 Pynamic Structure Factor, version 2.0

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






