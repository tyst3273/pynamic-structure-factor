import numpy as np
import os
import sys
sys.path.append('modules')
import Plot


Q = '3'

high_cut = 30 # meV
low_cut = 3 # meV
sigma = 0 # units = ???

f_dir = 'sqe/300K/constQ/'





##########################################################
################# avg and plot below here ################
##########################################################

### the Q_point for the plot label
Q_i = {'0':[2,5,5],
       '1':[5,3,5],
       '2':[6,0,1],
       '3':[3,2,3],
       '4':[7,0,0],
       '5':[7,1,0],
       '6':[2,3,4],
       '7':[2,4,3]}

### the file index to find the data in
f_i = {'0':'1',
       '1':'1',
       '2':'1',
       '3':'2',
       '4':'2',
       '5':'2',
       '6':'3',
       '7':'3'}

### position in file
p_i = {'0':0,
       '1':1,
       '2':2,
       '3':0,
       '4':1,
       '5':2,
       '6':0,
       '7':1}

nQ = 51 # number of Q per zone center

### get the right files
f_list = list(os.listdir(f_dir))
f_names = []
for f in f_list:
    if f.strip().split('_')[3] == f_i[Q]:
        f_names.append(f)

### get the Q_points from the Q_points file
inds = [0]
Q_points = np.loadtxt('Qpoints_'+f_i[Q])
Q_points = Q_points[nQ*p_i[Q]:nQ*(p_i[Q]+1),:]
for q in range(Q_points.shape[0]):
    res = Q_points[q,0]%1+Q_points[q,1]%1+Q_points[q,2]%1
    if res != 0:
        inds.append(q+1+nQ*p_i[Q])

### load and average the data
sqw = np.loadtxt(f_dir+f_names[0])
energy = sqw[:,0]
num_e = energy.shape[0]
if len(f_names) > 1:
    for f in f_names[1:]:
        sqw = sqw + np.loadtxt(f_dir+f)
sqw = sqw/len(f_names)
sqw = sqw[:,inds] # get the Q_points from above
sqw = sqw[:,1:].sum(axis=1)
sqw = np.append(energy.reshape((num_e,1)),sqw.reshape((num_e,1)),axis=1)

### trim the energies to the specified window
sqw = sqw[:sqw.shape[0]//2,:]
ind = np.argwhere(sqw[:,0] <= high_cut)[-1,0]
sqw = sqw[:ind,:]
ind = np.argwhere(sqw[:,0] >= low_cut)[0,0]
sqw = sqw[ind:,:]

### plot it
plot = Plot.plot_single_from_txt(sqw,Q_i[Q],sigma=sigma,fig_name=False,log_scale=False)





