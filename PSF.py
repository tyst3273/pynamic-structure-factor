from timeit import default_timer as timer
from mpi4py import MPI
import numpy as np
import sys

mod_path = '/home/ty/research/repos/pynamic-structure-factor/modules/'
sys.path.append(mod_path)

import mod_invars
import mod_io
import mod_lattice
import mod_Qpoints
import mod_sqw
import mod_utils 



# get the input file from command line
if len(sys.argv) != 1:
    input_file = sys.argv[1]
else:
    input_file = 'input_params'



# initialize the mpi stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
num_ranks = comm.Get_size()



if rank == 0:
    
    # start a timer
    start_time = timer()

    # print herald
    mod_utils.print_herald(num_ranks)

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

else:

    # receive the vars from rank 0
    invars = comm.recv(source=0,tag=0)
    lattice = comm.recv(source=0,tag=1)
    Qpoints = comm.recv(source=0,tag=2)



# get Qpoints for this rank, conver to 1/A
Qpoints.rank_init(lattice,rank)

# open the hdf5 file
traj_file = mod_io.traj_file(invars)

# setup calculator
sqw = mod_sqw.sqw(invars,Qpoints,rank)

# run the calculation
sqw.calculate(invars,Qpoints,lattice,traj_file)

# close the hdf5 files
traj_file.close()



if rank == 0:

    # gather all of the results
    sqw_from_ranks = [sqw.sqw]
    for ii in range(1,num_ranks):
        sqw_ii = comm.recv(source=ii,tag=11)
        sqw_from_ranks.append(sqw_ii)

    # assemble the results back into one 
    sqw_total = np.zeros((sqw.num_freq,Qpoints.total_Qsteps))
    mod_utils.assemble_sqw(sqw_from_ranks,sqw_total)

    # save it
    f_name = invars.outfile_prefix+f'_FINAL.hdf5'
    mod_io.save_sqw(invars,Qpoints.total_reduced_Q,sqw.meV,sqw_total,f_name)

    # calculate and print elapsed time
    end_time = timer()
    elapsed_time = (end_time-start_time)/60 #minutes
    message = f'total elapsed time: {elapsed_time:2.3f} minutes'
    mod_io.print_stdout(message,msg_type='TIMING')

else:

    # send the results to rank 0
    comm.send(sqw.sqw,dest=0,tag=11)



