
import numpy as np
import matplotlib.pyplot as plt

import psf.m_io as m_io



shape = [200,120]
sq = np.zeros(shape)

n_reps = 11
for ii in range(n_reps):
    
    f_name = f'./out/o{ii:g}_STRUFACS.hdf5'

    reader = m_io.c_reader(f_name)
    reader.read_diffuse()

    H = reader.H; K = reader.K; L = reader.L
    sq += reader.sq_diffuse[:,0,:]

sq /= n_reps


vmin = 0
vmax = 2.5e5
cmap = 'magma'
interp = 'none'


fig, ax = plt.subplots(figsize=(7,6))

extent = [H.min(),H.max(),L.min(),L.max()]

im = ax.imshow(sq,origin='lower',aspect='auto',extent=extent,cmap=cmap,
        interpolation=interp,vmin=vmin,vmax=vmax)
fig.colorbar(im,ax=ax,extend='both')

ax.axis(extent)
plt.show()



