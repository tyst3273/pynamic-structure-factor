import sys
import numpy as np
from mpi4py import MPI

sys.path.append('modules')

import Parser
import Parameters 
import FileIO
import Calculator
import ParalUtils 


# ==============================================================
# we initialize the communicator and get some data about the MPI
# options. rank is which process. 0 is the primary rank. 
# n_ranks is how many processes we are usign. 
# ==============================================================

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_ranks = comm.Get_size()


# ==============================================================
# rank 0 parses the input file, initializes stuff common to all
# ranks, and passes that around. rank_0_init generates the total
# list of Q points that will be split among processes. note, we 
# have to broadcast params before calling rank_x_init. this 
# method opens the hdf5 file and for some reason, we cant pass
# open files between processes.
# =============================================================

if rank == 0:

    FileIO.print_herald(n_ranks)

    parser = Parser.parser()
    parser.parse('input_params')

    params = Parameters.params(parser)
    params.rank_0_init()

    total_Qsteps = params.total_Qsteps
    total_reduced_Q = params.total_reduced_Q

    # will eventuall pass seperate chunks of Q points to be 
    # simultaneously done on each process. for now, this is a 
    # place holder.
    reduced_Q = total_reduced_Q 

    for ii in range(1,n_ranks):
        comm.send(reduced_Q,dest=ii,tag=0)
        comm.send(params,dest=ii,tag=1)

#    params.rank_x_init(reduced_Q,rank)
    params.rank_x_init(np.flipud(reduced_Q),rank)


# ====================================================================
# the other ranks receive the params object then each call rank_x_init 
# to do the rest of the initialization. once this is done for all 
# ranks, we can go to the calculator.
# ====================================================================

else:
   
    reduced_Q = comm.recv(source=0,tag=0)
    params = comm.recv(source=0,tag=1)

    params.rank_x_init(reduced_Q,rank)



# ===================================================================
# can now run the calculator on each process. in principle, this will
# be on a different set of Q-points on each process. for now, i am
# testing it on the same set of Q-points on all processes
# ===================================================================

calc = Calculator.calc(params)
calc.run(params)

# close the hdf5 file
params.clean_up()


# ===================================================================
# the calculator saves the data after each block and after all blocks.
# they are labeled *_P{rank}_B{block}.dat. i will figure out a way 
# to gather the sqw from each proc back into a single array and save 
# that one, but for now, each proc saves its own data. 
# ===================================================================



    
    

