import numpy as np
import os
import sys
sys.path.append('modules')
import Plot


high_cut = 50 # meV 
low_cut = 0 # meV
sigma = 0 # units = ???


f_dir = 'sqe/300K/constQ/'
f_names = list(os.listdir(f_dir))


sqw = np.loadtxt(f_dir+f_names[0])
energy = sqw[:,0]
num_e = energy.shape[0]
if len(f_names) > 1:
    for f in f_names[1:]:
        sqw = sqw + np.loadtxt(f_dir+f)
sqw = sqw/len(f_names)
sqw = sqw[:,1:].sum(axis=1)
sqw = np.append(energy.reshape((num_e,1)),sqw.reshape((num_e,1)),axis=1)


sqw = sqw[:sqw.shape[0]//2,:]
ind = np.argwhere(sqw[:,0] <= high_cut)[-1,0]
sqw = sqw[:ind,:]
ind = np.argwhere(sqw[:,0] >= low_cut)[0,0]
sqw = sqw[ind:,:]


plot = Plot.plot_single_from_txt(sqw,[2,3,4],sigma=sigma,fig_name='234_binned',log_scale=False)





