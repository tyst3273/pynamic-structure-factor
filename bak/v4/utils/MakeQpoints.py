"""
Input file is a csv with Qpoints in the format
Q1x Q1y Q1z
Q2x Q2y Q2z
...
"""

import numpy as np

in_file = 'Q_input'
out_file = 'Qpoints'


Q_input = np.loadtxt('Q_input')
if len(Q_input.shape) == 1:
    Q_input = Q_input.reshape((1,3))
num_Q = Q_input.shape[0]

Q_points = np.zeros((1,3))



dh = 0.05
dk = 0.05
dl = 0.1

for d0 in [0*dh,1*dh,-1*dh]:
    for d1 in [0*dk,1*dk,-1*dk]:
        for d2 in [0*dl,1*dl,-1*dl]:
            tmp_Q = np.copy(Q_input)
            tmp_Q[:,0] = tmp_Q[:,0]+d0
            tmp_Q[:,1] = tmp_Q[:,1]+d1
            tmp_Q[:,2] = tmp_Q[:,2]+d2
            Q_points = np.append(Q_points,tmp_Q,axis=0)




dh = 0.1
dk = 0.1
dl = 0.1

for d0 in [0*dh,1*dh,-1*dh]:
    for d1 in [0*dk,1*dk,-1*dk]:
        for d2 in [0*dl,1*dl,-1*dl]:
            tmp_Q = np.copy(Q_input)
            tmp_Q[:,0] = tmp_Q[:,0]+d0
            tmp_Q[:,1] = tmp_Q[:,1]+d1
            tmp_Q[:,2] = tmp_Q[:,2]+d2
            Q_points = np.append(Q_points,tmp_Q,axis=0)



Q_points = Q_points[1:,:]
Q_list = []
for Q in range(Q_points.shape[0]):
    Q_list.append(f'{Q_points[Q,0]:3.4f} {Q_points[Q,1]:3.4f} {Q_points[Q,2]:3.4f}')
Q_list = np.unique(np.array(Q_list))


np.savetxt(out_file,Q_list,fmt='%s')








