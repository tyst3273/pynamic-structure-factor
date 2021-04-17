import sys
import numpy as np
from mpi4py import MPI

sys.path.append('modules')

import Parser
import Parameters 
import FileIO
import Calculator
from ParaUtils import *


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_ranks = comm.Get_size()


if rank == 0:

    # print the herald 
    FileIO.print_herald(n_ranks)

    # get input opts from file
    parser = Parser.parser()
    parser.parse('input_params')

    # initialize params obj
    params = Parameters.params(parser)
    params.rank_0_init()

    # get the total set of Q's 
    total_Qsteps = params.total_Qsteps
    total_reduced_Q = params.total_reduced_Q

    # split the total set of Q's over the procs
    prepare_Qpoints(total_reduced_Q,total_Qsteps,n_ranks,comm)

    reduced_Q = total_reduced_Q 

    # pass the sets of Q's and params object to each proc.
    for ii in range(1,n_ranks):
        comm.send(reduced_Q,dest=ii,tag=0)
        comm.send(params,dest=ii,tag=1)

    # initiliaze params on this proc using set of Q's
    params.rank_x_init(np.flipud(reduced_Q),rank)


else:
   
    # get the set of Q's and params object
    reduced_Q = comm.recv(source=0,tag=0)
    params = comm.recv(source=0,tag=1)

    # initiliaze params on this proc using set of Q's
    params.rank_x_init(reduced_Q,rank)


comm.Abort()

# run the calculation on all procs
calc = Calculator.calc(params)
calc.run(params)

# close the hdf5 files
params.clean_up()



    
    

