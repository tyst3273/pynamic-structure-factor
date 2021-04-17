import numpy as np
from FileIO import print_stdout

def prepare_Qpoints(total_reduced_Q,total_Qsteps,n_ranks,comm):

    if n_ranks == 1:
        return [total_reduced_Q]

    if total_Qsteps < n_ranks:
        if n_ranks == 1:
            message = f'cannot split {total_Qsteps} over {n_ranks}. aborting'
            print_stdout(message,msg_type='ERROR')
            exit()
        else:
            message = f'cannot split {total_Qsteps} Q-points over {n_ranks} processors. aborting'
            print_stdout(message,msg_type='ERROR')
            comm.Abort()

    num_Q_per_proc = np.zeros(n_ranks)

    npp = total_Qsteps // n_ranks
    npp_remainder = total_Qsteps % n_ranks
    
    num_Q_per_proc[:] = npp
    num_Q_per_proc[0] = num_Q_per_proc[0]+npp_remainder

    if npp_remainder != 0:
        message = f' one proc. will compute {npp_remainder+npp} Qpoints; the rest will compute {npp}'
        print_stdout(message,msg_type='PARALLELISM')

    else:
        message = f' each proc. will compute {npp} Qpoints'
        print_stdout(message,msg_type='PARALLELISM')

    print(num_Q_per_proc)






