
import matplotlib.pyplot as plt
import numpy as np

import psf.m_io as m_io


T = 300
kB = 0.08617333262

reader = m_io.c_reader('cmap_STRUFACS.hdf5')
reader.read_sqw()

Q = reader.Q_rlu
nQ = Q.shape[0]
sqw = np.fft.fftshift(reader.sqw,axes=(1))
energy = np.fft.fftshift(reader.energy)
extent = [0,1,energy.min(),energy.max()]

fig, ax = plt.subplots(figsize=(12,6))

#bose = np.exp(energy/(kB*T))-1

im = ax.imshow(sqw.T,aspect='auto',origin='lower',cmap='hot',extent=extent,
    vmin=0,vmax=50)
fig.colorbar(im,ax=ax,extend='both')

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='medium')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
ax.axis([0,1,0,75])

ax.set_xlabel(r'Q',labelpad=1.0,fontsize='large')
ax.set_ylabel(r'Energy [meV]',labelpad=1.0,fontsize='large')

plt.savefig('lammps_sqw.pdf',bbox_inches='tight')


