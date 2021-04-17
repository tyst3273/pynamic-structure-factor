import numpy as np
import os
import sys
sys.path.append('modules')
import Plot


high_cut = 40 # meV 
low_cut = 10 # meV
sigma = 2 # units = ???


f_dir = 'sqe/300K/constQ/'
f_names = list(os.listdir(f_dir))


sqw = np.loadtxt(f_dir+f_names[0])
energy = sqw[:,0]
if len(f_names) > 1:
    for f in f_names[1:]:
        sqw = sqw + np.loadtxt(f_dir+f)
sqw = sqw/len(f_names)
sqw[:,0] = energy


sqw = sqw[:sqw.shape[0]//2,:]
ind = np.argwhere(sqw[:,0] <= high_cut)[-1,0]
sqw = sqw[:ind,:]
ind = np.argwhere(sqw[:,0] >= low_cut)[0,0]
sqw = sqw[ind:,:]


with open(f_dir+f_names[0],'r') as fid:
    Qpoints = np.array(fid.readline().strip().strip('#').split()[1:]).astype(float)
Qpoints = Qpoints.reshape((Qpoints.shape[0]//3,3))

for Q in range(Qpoints.shape[0]):
    Qpoint = Qpoints[Q,:]
    plot = Plot.plot_single_from_txt(sqw[:,[0,Q+1]],Qpoint,sigma=sigma,fig_name=False,log_scale=False)





