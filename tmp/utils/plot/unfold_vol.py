import os
import matplotlib.pyplot as plt
import numpy as np
import sys
import h5py 
from timeit import default_timer
from scipy.interpolate import interpn
from libpsf.mod_io import read_sqw, read_bragg, read_timeavg

# --------------------------------------------------------------------------------------------------
class timer:

    def __init__(self,label,units='s'):
    
        if units == 'ms':
            self.scale = 1e3
            self.units = 'ms'
        elif units == 'm':
            self.scale = 1/60
            self.units = 'm'
        else:
            self.scale = 1
            self.units = 's'

        self.label = label

    def start(self):
        self.start_time = default_timer()

    def stop(self):
        self.stop_time = default_timer()
        self.lap_time = (self.stop_time-self.start_time)*self.scale
        print(f'\n timing info from timer \'{self.label}\'\n   lap time: {self.lap_time:<8.4f}'\
              f' [{self.units}]\n')


# --------------------------------------------------------------------------------------------------

# Q range in the data file
Q_range = [0,3]
dQ = 0.1
num_Q = 31

# xray
energy, Qpts, raw = read_sqw('ins_333/deut_300_xray_333_SQW_FINAL.hdf5')
Qpts = np.array(Qpts,dtype=float)
Qpts = np.round(Qpts,1)

print(Qpts)
# get the inds that fall inside the energy range
e_int = 1
e_inds = np.argwhere(energy >= -e_int).flatten()
e_inds = np.intersect1d(np.argwhere(energy <= e_int).flatten(),e_inds)

# slice the data
raw = raw[e_inds,:]
raw = np.mean(raw,axis=0)

# Q array
Q_arr = np.arange(Q_range[0],Q_range[1]+dQ,dQ)
full_Q_arr = np.arange(-1*Q_range[1],Q_range[1]+dQ,dQ)
Q_arr = np.round(Q_arr,2)
full_Q_arr = np.round(full_Q_arr,2)

# unfolded onto positive quadrant of BZ grid 
pos_bz = np.zeros((num_Q,num_Q,num_Q)) 
num_Q_file = Qpts.shape[0]

loop_timer = timer('loop_timer')
loop_timer.start()

pos_bz = np.zeros((num_Q,num_Q,num_Q))
for qq in range(num_Q_file):
    Q = Qpts[qq,:]
    hf = np.modf(Q[0])[0]
    kf = np.modf(Q[1])[0]
    lf = np.modf(Q[2])[0]
    h_ind = np.argwhere(Q_arr == Q[0]).flatten()[0]
    k_ind = np.argwhere(Q_arr == Q[1]).flatten()[0]
    l_ind = np.argwhere(Q_arr == Q[2]).flatten()[0]
    inds = np.zeros((6,3),dtype=int)
    inds[0,:] = [h_ind,k_ind,l_ind]
    inds[1,:] = [h_ind,l_ind,k_ind]
    inds[2,:] = [k_ind,h_ind,l_ind]
    inds[3,:] = [k_ind,l_ind,h_ind]
    inds[4,:] = [l_ind,h_ind,k_ind]
    inds[5,:] = [l_ind,k_ind,h_ind]
    for ii in range(6):
        pos_bz[inds[ii,0],inds[ii,1],inds[ii,2]] = raw[qq]

print(pos_bz.shape)
# write it to an hdf5 file
with h5py.File('md_xray_pos_octant.hdf5','w') as out_db:
    Q = out_db.create_dataset('Q_arr',[num_Q])
    Q[:] = Q_arr[:]
    I = out_db.create_dataset('intensity',[num_Q,num_Q,num_Q])
    I[:,:,:] = pos_bz[:,:,:]
loop_timer.stop()

# interpolate into fine grid
q_lo = 0; q_hi = 3
n_fine = 81
q = np.linspace(q_lo,q_hi,n_fine)
x, y, z = np.meshgrid(q,q,q)
nq = x.size
x.shape = [nq,1]
y.shape = [nq,1]
z.shape = [nq,1]
points = np.append(x,y,axis=1)
points = np.append(points,z,axis=1)

fine_grid = interpn((Q_arr,Q_arr,Q_arr),pos_bz,points)
fine_grid.shape = [n_fine,n_fine,n_fine]

# write it to an hdf5 file
with h5py.File('md_xray_pos_bz_fine.hdf5','w') as out_db:
    Q = out_db.create_dataset('Q_arr',[n_fine])
    Q[:] = q[:]
    I = out_db.create_dataset('intensity',[n_fine,n_fine,n_fine])
    I[:,:,:] = fine_grid[:,:,:]


# -------------------------------------------------------------------------

# now unfold onto full bz
full_bz = np.zeros((2*num_Q-1,2*num_Q-1,2*num_Q-1))
full_bz[:num_Q,:num_Q,:num_Q] = np.flip(pos_bz,axis=[0,1,2])
full_bz[:num_Q,num_Q-1:,:num_Q] = np.flip(pos_bz,axis=[0,2])
full_bz[:num_Q,:num_Q,num_Q-1:] = np.flip(pos_bz,axis=[0,1])
full_bz[:num_Q,num_Q-1:,num_Q-1:] = np.flip(pos_bz,axis=[0])
full_bz[num_Q-1:,:num_Q,:num_Q] = np.flip(pos_bz,axis=[1,2])
full_bz[num_Q-1:,num_Q-1:,:num_Q] = np.flip(pos_bz,axis=[2])
full_bz[num_Q-1:,:num_Q,num_Q-1:] = np.flip(pos_bz,axis=[1])
full_bz[num_Q-1:,num_Q-1:,num_Q-1:] = pos_bz
full_Q = np.zeros(2*num_Q-1)
full_Q[:num_Q] = -np.flip(Q_arr)
full_Q[num_Q:] = Q_arr[1:]

if False:
    # the l to plot
    which_l = 0.5
    # get the index
    which_l = np.round(float(which_l),2)
    if which_l in full_Q_arr:
        l_ind = np.argwhere(full_Q_arr == which_l).flatten()[0]
    else:
        print('fuck!')
    # plot the bz slice
    fig, ax = plt.subplots(figsize=[6,6])
    im = ax.imshow(full_bz[:,:,l_ind],origin='lower',aspect='auto',vmax=5e3,interpolation='nearest')
    fig.colorbar(im,ax=ax)
    plt.show()

num_Q = full_bz.shape[0]
# write it to an hdf5 file
with h5py.File('md_xray_full_bz.hdf5','w') as out_db:
    Q = out_db.create_dataset('Q_arr',[num_Q])
    Q[:] = full_Q[:]
    I = out_db.create_dataset('intensity',[num_Q,num_Q,num_Q])
    I[:,:,:] = full_bz[:,:,:]

# interpolate into fine grid
q_lo = -2; q_hi = 2
n_fine = 101
q = np.linspace(q_lo,q_hi,n_fine)
x, y, z = np.meshgrid(q,q,q)
nq = x.size
x.shape = [nq,1]
y.shape = [nq,1]
z.shape = [nq,1]
points = np.append(x,y,axis=1)
points = np.append(points,z,axis=1)

fine_grid = interpn((full_Q,full_Q,full_Q),full_bz,points)
fine_grid.shape = [n_fine,n_fine,n_fine]

if False:
    # the l to plot
    which_l = 0.5
    # get the index
    which_l = np.round(float(which_l),1)
    if which_l in q:
        l_ind = np.argwhere(q == which_l).flatten()[0]
    else:
        print('fuck!')
    # plot the bz slice
    fig, ax = plt.subplots(figsize=[6,6])
    im = ax.imshow(fine_grid[:,:,l_ind],origin='lower',aspect='auto',vmax=5e3,interpolation='none')
    fig.colorbar(im,ax=ax)
    plt.show()


num_Q = n_fine
# write it to an hdf5 file
with h5py.File('md_xray_full_bz_fine.hdf5','w') as out_db:
    Q = out_db.create_dataset('Q_arr',[num_Q])
    Q[:] = q[:]
    I = out_db.create_dataset('intensity',[num_Q,num_Q,num_Q])
    I[:,:,:] = fine_grid[:,:,:]


















