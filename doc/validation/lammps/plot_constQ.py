
import matplotlib.pyplot as plt
import numpy as np

import psf.m_io as m_io


T = 300
kB = 0.08617333262

reader = m_io.c_reader('constQ_STRUFACS.hdf5')
reader.read_sqw()

Q = reader.Q_rlu
nQ = Q.shape[0]
sqw = np.fft.fftshift(reader.sqw,axes=(1))
energy = np.fft.fftshift(reader.energy)

fig, ax = plt.subplots(figsize=(9,9))

#bose = np.exp(energy/(kB*T))-1

d = 10
s = 0
for qq in range(nQ):
    ax.plot(energy,sqw[qq,:]+s,marker='o',ms=0,
        lw=2,label=f'Q=({Q[qq,0]:3.2f},{Q[qq,1]:3.2f},{Q[qq,2]:3.2f})')
    s += d

fig.legend()


for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='medium')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
ax.axis([-75,75,-5,100])

ax.set_ylabel(r'Intensity [arb. units]',labelpad=1.0,fontsize='large')
ax.set_xlabel(r'Energy [meV]',labelpad=1.0,fontsize='large')

plt.savefig('lammps_constQ.pdf',bbox_inches='tight')


