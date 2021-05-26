# creats a set of Q points for 'binning' around a Q vector in reciprocal space.
# this isn't complicated, adapt as needed or write your own tool.

import numpy as np


### in and out files
in_file = 'input_Qs'
out_file = 'Qi_'


### loop over input Q's and append to the output list
Q_input = np.loadtxt(in_file)
num_Q = Q_input.shape[0]



for Q in range(num_Q):

    Q_output = np.zeros((1,3))

    Qi = Q_input[Q,:]
    print(f'Q = {Qi[0]:2.2f} {Qi[1]:2.2f} {Qi[2]:2.2f}')
    
    ### bin by 0.05
    dh = 0.05
    dk = 0.05
    dl = 0.1

    for d0 in [0*dh,1*dh,-1*dh]:
        for d1 in [0*dk,1*dk,-1*dk]:
            for d2 in [0*dl,1*dl,-1*dl]:
                tmp_Q = np.copy(Qi).reshape((1,3))
                tmp_Q[:,0] = tmp_Q[:,0]+d0
                tmp_Q[:,1] = tmp_Q[:,1]+d1
                tmp_Q[:,2] = tmp_Q[:,2]+d2
                Q_output = np.append(Q_output,tmp_Q,axis=0)

    ### bin by 0.1
    dh = 0.1
    dk = 0.1
    dl = 0.1

    for d0 in [0*dh,1*dh,-1*dh]:
        for d1 in [0*dk,1*dk,-1*dk]:
            for d2 in [0*dl,1*dl,-1*dl]:
                tmp_Q = np.copy(Qi).reshape((1,3))
                tmp_Q[:,0] = tmp_Q[:,0]+d0
                tmp_Q[:,1] = tmp_Q[:,1]+d1
                tmp_Q[:,2] = tmp_Q[:,2]+d2
                Q_output = np.append(Q_output,tmp_Q,axis=0)

    ### trim the 0,0,0 from the output list
    Q_output = Q_output[1:,:]
    Q_list = []

    ### format as string to make unique
    for i in range(Q_output.shape[0]):
        Q_list.append(f'{Q_output[i,0]:3.4f} {Q_output[i,1]:3.4f} {Q_output[i,2]:3.4f}')
    Q_list = np.unique(np.array(Q_list))

    np.savetxt(out_file+f'{Q}',Q_list,fmt='%s')






