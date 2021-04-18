import sys
import numpy as np
from timeit import default_timer as timer
from mpi4py import MPI

sys.path.append('/home/ty/custom_modules/pynamic-structure-factor/')

import Parser
import Parameters 
import FileIO
import Calculator
from ParalUtils import *


# initialize the communicator and get info
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_ranks = comm.Get_size()


if rank == 0:

    # start a timer
    start_time = timer()

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
    reduced_Q_set = prepare_Qpoints(params,n_ranks)

    # pass the sets of Q's and params object to each proc.
    for ii in range(1,n_ranks):
        comm.send(reduced_Q_set[ii],dest=ii,tag=0)
        comm.send(params,dest=ii,tag=1)

    # initiliaze params on rank 0 using set of Q's
    params.rank_x_init(reduced_Q_set[0],rank)

else:
   
    # get the set of Q's and params object
    reduced_Q = comm.recv(source=0,tag=0)
    params = comm.recv(source=0,tag=1)

    # initiliaze params on this proc using set of Q's
    params.rank_x_init(reduced_Q,rank)


# initialize calculation on all procs
calc = Calculator.calc(params)

# run calc on all procs and close the hdf5 file when done. 
calc.run(params)


if rank == 0:
    
    # gather all of the results
    sqw_total_set = [calc.sqw]
    for ii in range(1,n_ranks):
        sqw_total_set.append(comm.recv(source=ii,tag=3))
    
    # assemble the results back into one
    sqw_total = assemble_sqw_total(params,sqw_total_set,n_ranks)
    
    # save it
    f_name = params.outfile_prefix+f'_FINAL.dat'
    FileIO.save_sqw(params,params.total_reduced_Q,sqw_total,f_name=f_name)

    # calculate and print elapsed time
    end_time = timer()
    elapsed_time = (end_time-start_time)/60 #minutes
    message = f'total elapsed time: {elapsed_time:2.3f} minutes'
    FileIO.print_stdout(message,msg_type='TIMING')
 
else:    

    # send the results to rank 0
    comm.send(calc.sqw,dest=0,tag=3)
    








