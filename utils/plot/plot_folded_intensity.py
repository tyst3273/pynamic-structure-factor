import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.io import loadmat
from matplotlib.image import imread
import sys
sys.path.append('/home/ty/research/repos/pynamic-structure-factor/modules')
from mod_io import read_sqw, read_bragg, read_timeavg

# -----------------------------------------------------------------------

# cut off range range of averaging SQW
e_cut = 2 # set to 2 to remake example low_energy_integrated.pdf

# options for plotting
num_Q = 41 # number of steps along each direction in the positive quadrant
Q_range = [0,2, 0,2] # range of Q in the positive quadrant

extent = [-2,2,-2,2] # extent of unfolded Q

interp = 'none' # colormap interpolation
cmap = 'CMRmap' # colormap

# ---------------------------------------------------------------------

num_Q = 41 #number of steps along path in the file
Q_range = [0,2, 0,2] # in the file

extent = [-2,2,-2,2]
vlims = [0,0.015]
#vlims = [None,None]

interp = 'none'
cmap = 'CMRmap'

# load sqw_raw
in_file = f'sqw/diffuse_BRAGG_FINAL.hdf5'
Qpts, bragg_raw = read_bragg(in_file)
in_file = f'sqw/diffuse_TIMEAVG_FINAL.hdf5'
Qpts, timeavg_raw = read_timeavg(in_file)
in_file = f'sqw/diffuse_SQW_FINAL.hdf5'
energy, Qpts, sqw_raw = read_sqw(in_file)

e_cut = 2 # meV
e_inds = np.argwhere(energy <= e_cut).flatten()
e_inds = np.intersect1d(e_inds,np.argwhere(energy >= -e_cut).flatten())
sqw_low_e_raw = sqw_raw[e_inds,:]
sqw_low_e_raw = sqw_low_e_raw.mean(axis=0)

e_cut = 100 # the max in dataset is ~|42| meV
e_inds = np.argwhere(energy <= e_cut).flatten()
e_inds = np.intersect1d(e_inds,np.argwhere(energy >= -e_cut).flatten())
sqw_raw = sqw_raw[e_inds,:]
sqw_raw = sqw_raw.mean(axis=0)

# ---------------------------------------------------------------------

# fill folded arrays
h_arr = np.linspace(Q_range[0],Q_range[1],num_Q)
k_arr = np.linspace(Q_range[2],Q_range[3],num_Q)

bragg_folded = np.zeros((num_Q,num_Q))
timeavg_folded = np.zeros((num_Q,num_Q))
sqw_folded = np.zeros((num_Q,num_Q))
sqw_low_e_folded = np.zeros((num_Q,num_Q))
for Qh in range(num_Q):
    for Qk in range(num_Q):
        h = np.round(h_arr[Qh],3) # precision of the Q points in file
        k = np.round(k_arr[Qk],3) 
        hi = np.argwhere(Qpts[:,0]==h).flatten()
        ki = np.argwhere(Qpts[:,1]==k).flatten()
        ind = np.intersect1d(hi,ki).flatten()
        if ind.shape[0] == 0:
            continue
        else:
            timeavg_folded[Qh,Qk] = timeavg_raw[ind]
            bragg_folded[Qh,Qk] = bragg_raw[ind]
            sqw_folded[Qh,Qk] = sqw_raw[ind]
            sqw_low_e_folded[Qh,Qk] = sqw_low_e_raw[ind]

# ----------------------------------------------------------------------

# unfold time averaged intensity onto full x-y plane
tmp = np.copy(timeavg_folded)
np.fill_diagonal(tmp,0)
timeavg_folded = timeavg_folded+tmp.T

timeavg = np.zeros((num_Q*2-1,num_Q*2-1))
timeavg[:num_Q,:num_Q] = np.fliplr(np.flipud(timeavg_folded[:,:]))
timeavg[:num_Q,num_Q:] = np.flipud(timeavg_folded[:,1:])
timeavg[num_Q:,:num_Q] = np.fliplr(timeavg_folded[1:,:])
timeavg[num_Q:,num_Q:] = timeavg_folded[1:,1:]

# unfold bragg intensity onto full x-y plane
tmp = np.copy(bragg_folded)
np.fill_diagonal(tmp,0)
bragg_folded = bragg_folded+tmp.T

bragg = np.zeros((num_Q*2-1,num_Q*2-1))
bragg[:num_Q,:num_Q] = np.fliplr(np.flipud(bragg_folded[:,:]))
bragg[:num_Q,num_Q:] = np.flipud(bragg_folded[:,1:])
bragg[num_Q:,:num_Q] = np.fliplr(bragg_folded[1:,:])
bragg[num_Q:,num_Q:] = bragg_folded[1:,1:]

# unfold sqw intensity onto full x-y plane
tmp = np.copy(sqw_folded)
np.fill_diagonal(tmp,0)
sqw_folded = sqw_folded+tmp.T

sqw = np.zeros((num_Q*2-1,num_Q*2-1))
sqw[:num_Q,:num_Q] = np.fliplr(np.flipud(sqw_folded[:,:]))
sqw[:num_Q,num_Q:] = np.flipud(sqw_folded[:,1:])
sqw[num_Q:,:num_Q] = np.fliplr(sqw_folded[1:,:])
sqw[num_Q:,num_Q:] = sqw_folded[1:,1:]

tmp = np.copy(sqw_low_e_folded)
np.fill_diagonal(tmp,0)
sqw_low_e_folded = sqw_low_e_folded+tmp.T

sqw_low_e = np.zeros((num_Q*2-1,num_Q*2-1))
sqw_low_e[:num_Q,:num_Q] = np.fliplr(np.flipud(sqw_low_e_folded[:,:]))
sqw_low_e[:num_Q,num_Q:] = np.flipud(sqw_low_e_folded[:,1:])
sqw_low_e[num_Q:,:num_Q] = np.fliplr(sqw_low_e_folded[1:,:])
sqw_low_e[num_Q:,num_Q:] = sqw_low_e_folded[1:,1:]

# -------------------------------------------------------------------------

vlims = [0,0.001] # colormap scale
fig,ax=plt.subplots(figsize=(7,6))
im = ax.imshow(timeavg-bragg,aspect='auto',cmap=cmap,origin='lower',extent=extent,
        interpolation=interp,vmin=vlims[0],vmax=vlims[1])
fig.colorbar(im,ax=ax,extend='both')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
xlabel = r'$\xi$ (r.l.u.)'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = r'$\eta$ (r.l.u.)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')
fig.suptitle(r"$\bf{Q}$=($\xi$,$\eta$,2)",fontsize='x-large',y=0.95)
plt.savefig(f'diffuse.pdf',format='pdf',dpi=150,bbox_inches='tight')
#plt.show()

vlims = [0,0.001] # colormap scale
fig,ax=plt.subplots(figsize=(7,6))
im = ax.imshow(bragg,aspect='auto',cmap=cmap,origin='lower',extent=extent,
        interpolation=interp,vmin=vlims[0],vmax=vlims[1])
fig.colorbar(im,ax=ax,extend='both')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
xlabel = r'$\xi$ (r.l.u.)'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = r'$\eta$ (r.l.u.)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')
fig.suptitle(r"$\bf{Q}$=($\xi$,$\eta$,2)",fontsize='x-large',y=0.95)
plt.savefig(f'bragg.pdf',format='pdf',dpi=150,bbox_inches='tight')
#plt.show()

vlims = [0,0.001] # colormap scale
fig,ax=plt.subplots(figsize=(7,6))
im = ax.imshow(timeavg,aspect='auto',cmap=cmap,origin='lower',extent=extent,
        interpolation=interp,vmin=vlims[0],vmax=vlims[1])
fig.colorbar(im,ax=ax,extend='both')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
xlabel = r'$\xi$ (r.l.u.)'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = r'$\eta$ (r.l.u.)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')
fig.suptitle(r"$\bf{Q}$=($\xi$,$\eta$,2)",fontsize='x-large',y=0.95)
plt.savefig(f'timeavg.pdf',format='pdf',dpi=150,bbox_inches='tight')
#plt.show()

vlims = [0,0.01] # colormap scale
fig,ax=plt.subplots(figsize=(7,6))
im = ax.imshow(sqw_low_e,aspect='auto',cmap=cmap,origin='lower',extent=extent,
        interpolation=interp,vmin=vlims[0],vmax=vlims[1])
fig.colorbar(im,ax=ax,extend='both')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
xlabel = r'$\xi$ (r.l.u.)'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = r'$\eta$ (r.l.u.)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')
fig.suptitle(r"$\bf{Q}$=($\xi$,$\eta$,2)",fontsize='x-large',y=0.95)
plt.savefig(f'low_energy_integrated.pdf',format='pdf',dpi=150,bbox_inches='tight')
#plt.show()

vlims = [0,0.001] # colormap scale
fig,ax=plt.subplots(figsize=(7,6))
im = ax.imshow(sqw,aspect='auto',cmap=cmap,origin='lower',extent=extent,
        interpolation=interp,vmin=vlims[0],vmax=vlims[1])
fig.colorbar(im,ax=ax,extend='both')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
xlabel = r'$\xi$ (r.l.u.)'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = r'$\eta$ (r.l.u.)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')
fig.suptitle(r"$\bf{Q}$=($\xi$,$\eta$,2)",fontsize='x-large',y=0.95)
plt.savefig(f'energy_integrated.pdf',format='pdf',dpi=150,bbox_inches='tight')
#plt.show()

vlims = [0,0.00001] # colormap scale
fig,ax=plt.subplots(figsize=(7,6))
im = ax.imshow(sqw-timeavg,aspect='auto',cmap=cmap,origin='lower',extent=extent,
        interpolation=interp,vmin=vlims[0],vmax=vlims[1])
fig.colorbar(im,ax=ax,extend='both')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
xlabel = r'$\xi$ (r.l.u.)'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = r'$\eta$ (r.l.u.)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')
fig.suptitle(r"$\bf{Q}$=($\xi$,$\eta$,2)",fontsize='x-large',y=0.95)
plt.savefig(f'energy_integrated_minus_timeavg.pdf',format='pdf',dpi=150,bbox_inches='tight')
#plt.show()


