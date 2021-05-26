
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
num_ranks = comm.Get_size()

