import numpy as np
import os
import sys
sys.path.append('modules')
import Plot


high_cut = 40 # meV 
low_cut = 0 # meV

f_dir = 'sqe/10K/2d/'
f_names = list(os.listdir(f_dir))

print(f_names[0])
sqw = np.loadtxt(f_dir+f_names[0])
energy = sqw[:,0]
if len(f_names) > 1:
    for f in f_names[1:]:
        print(f)
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
Qmax = Qpoints[-3:]
Qmin = Qpoints[:3]

plot = Plot.plot_from_txt(sqw,Qmin,Qmax,fig_name='test') #,vmin=-2,vmax=7,cmap='jet')





