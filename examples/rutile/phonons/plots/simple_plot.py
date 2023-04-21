
import numpy as np
import matplotlib.pyplot as plt

from psf.m_io import c_reader


reader = c_reader('./../out/pristine_exe_STRUFACS.hdf5')
reader.read_sqw()
Q_rlu = reader.Q_rlu
energy = reader.energy
sqw = reader.sqw

num_verts = 3
num_e = energy.size//2

cmap = 'viridis'
interp = 'none'
vmin = 0
vmax = 5e-5
c = (1,1,1)

fig, ax = plt.subplots(figsize=(7,4))


extent = [0,num_verts,0,energy.max()]

im = ax.imshow(sqw[:,:num_e].T,origin='lower',aspect='auto',interpolation=interp,cmap=cmap,
                 vmin=vmin,vmax=vmax,extent=extent)

fig.colorbar(im,ax=ax,location='right',extend='both',aspect=20,pad=0.03)


ax.set_ylim(0,120)
ax.set_xlim(0,3)
ax.set_xticks([0,1,2,3])

ax.plot([2,2],[0,energy.max()],lw=2,ls='-',c='k')
ax.plot([2,2],[0,energy.max()],lw=2,ls='-',c='k')

ax.set_xticklabels(['0,0,0','5,0,0','5,5,0; 2.5,0,0','2.5,5,0'])
ax.set_xlabel('Q (rlu)',labelpad=4,fontsize='large')


ax.set_ylabel('Energy (meV)',labelpad=4,fontsize='large')

plt.show()

