
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm 
from matplotlib.colors import Normalize
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"

from psf.m_io import c_reader


fig_name = 'test.png'
reader = c_reader('mesh_STRUFACS.hdf5')
reader.read_energy_integrated_sqw()
sq = reader.sq_integrated.squeeze()

H = reader.H
K = reader.K
L = reader.L

sq = np.fft.fftshift(sq,axes=1)
energy = reader.energy

cmap = 'viridis'
interp = 'none'
vmin = 0 #-1e-7
#vmax = 0.05
vmax = 0.1 #1e-7
c = (1,1,1)

fig, ax = plt.subplots(figsize=(6,5))

extent = [H.min(),H.max(),K.min(),K.max()]

norm = Normalize
im = ax.imshow(sq.T,origin='lower',aspect='auto',interpolation=interp,cmap=cmap,
                 extent=extent,norm=norm(vmin,vmax))
cbar = fig.colorbar(im,ax=ax,location='right',extend='both',aspect=40,pad=0.03)

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='medium')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
#ax.set_xticks([-3,-2,-1,0,1,2,3])
#ax.set_yticks([-3,-2,-1,0,1,2,3])
ax.set_rasterized(True)

ax.set_ylabel('H (rlu)',labelpad=4,fontsize='large')
ax.set_xlabel('K (rlu)',labelpad=4,fontsize='large')

plt.savefig(fig_name,dpi=150,bbox_inches='tight')

plt.show()

