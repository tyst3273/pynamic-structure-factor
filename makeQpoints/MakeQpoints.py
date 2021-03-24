"""
Input file is a csv with Qpoints in the format
Q1x Q1y Q1z
Q2x Q2y Q2z
...

Output is a csv file with the points around each Q grouped in order
of Q's in the unput file
"""

import numpy as np

### in and out files
in_file = 'Qi_1'
out_file = 'Qo_1'



### read the input data
Q_input = np.loadtxt(in_file)
if len(Q_input.shape) == 1:
    Q_input = Q_input.reshape((1,3))
num_Q = Q_input.shape[0]

### loop over input Q's and append to the output list
Q_output = np.zeros((1,3))

for Q in range(num_Q):
    
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
for Q in range(Q_output.shape[0]):
    Q_list.append(f'{Q_output[Q,0]:3.4f} {Q_output[Q,1]:3.4f} {Q_output[Q,2]:3.4f}')
Q_list = np.unique(np.array(Q_list))

np.savetxt(out_file,Q_list,fmt='%s')






