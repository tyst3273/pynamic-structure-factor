import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.io import loadmat
from matplotlib.image import imread
import sys

sys.path.append('/home/ty/custom_modules/pynamic-structure-factor/')

from mod_io import read_sqw


e_lo = 0
e_hi = 2

num_Q = 81
Q_range = [0,4,0,4]

extent = [-4,4,-4,4]
vlims = [0,15]
#vlims = [None,None]

interp = 'none'
cmap = 'hot'

which = 'all'
in_file = f'{which}_FINAL.hdf5'
in_file = os.path.join(os.getcwd(),'diffuse',in_file)


############ load data #######################

energy, qpts, dat = read_sqw(in_file)

num_e = energy.shape[0]//2
energy = energy[:num_e]
e_max = energy.max()
dat = dat[:num_e,:]

######### energy integral ########
e_lo_ind = np.argwhere(energy >= e_lo).flatten().min()
e_hi_ind = np.argwhere(energy <= e_hi).flatten().max()
e_slice = energy[e_lo_ind:e_hi_ind]
dat = dat[e_lo_ind:e_hi_ind,:]
dat = np.trapz(dat,e_slice,axis=0)

######## fill SQE array ##########
h_arr = np.linspace(Q_range[0],Q_range[1],num_Q)
k_arr = np.linspace(Q_range[2],Q_range[3],num_Q)

#qpts = np.loadtxt('qpt_L0.5')

sqe = np.zeros((num_Q,num_Q))
for Qh in range(num_Q):
    for Qk in range(num_Q):
        h = np.round(h_arr[Qh],3)
        k = np.round(k_arr[Qk],3)
        hi = np.argwhere(qpts[:,0]==h).flatten()
        ki = np.argwhere(qpts[:,1]==k).flatten()
        ind = np.intersect1d(hi,ki)
        if ind.shape[0] == 0:
            continue
        else:
            sqe[Qh,Qk] = dat[ind[0]]


tmp = np.copy(sqe)
np.fill_diagonal(tmp,0)
sqe = sqe+tmp.T


dat = np.zeros((num_Q*2,num_Q*2))
dat[:num_Q,:num_Q] = np.fliplr(np.flipud(sqe))
dat[:num_Q,num_Q:] = np.flipud(sqe)
dat[num_Q:,:num_Q] = np.fliplr(sqe)
dat[num_Q:,num_Q:] = sqe


######### create plots ###########
fig,ax=plt.subplots(figsize=(7,6))
im = ax.imshow(dat,aspect='auto',cmap=cmap,origin='lower',extent=extent,
        interpolation=interp,vmin=vlims[0],vmax=vlims[1])

fig.colorbar(im,ax=ax,extend='both')

######## format plots #########
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
#ax.axis([0,2,0,2])

xlabel = r'($\xi$,0,0) (r.l.u.)'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = r'(0,$\xi$,0) (r.l.u.)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')

ax.annotate(r"{} modes".format(which),xy=(0.02,0.97),
        xycoords="axes fraction",color='w',fontsize='large')

fig.suptitle(r"$\bf{Q}$=(H,K,0.5)",fontsize='x-large',y=0.95)

np.savetxt(f'HK0.5_{which}_integrated.dat',dat)
plt.savefig(f'HK0.5_{which}.pdf',format='pdf',dpi=150,bbox_inches='tight')

plt.show()



