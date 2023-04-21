
import numpy as np
import h5py
import matplotlib.pyplot as plt


f = 'exp_cuts.hdf5'

with h5py.File(f,'r') as db:

    x = db['xray_2.5'][...]

nQ = x.shape[0]

 
steps = 300

x = x[steps:nQ-steps,steps:nQ-steps]

plt.imshow(x,cmap='viridis',interpolation='none',aspect='auto',origin='lower',vmin=0,vmax=500)
plt.show()


np.savetxt('xds_cut_L2.5',x)



