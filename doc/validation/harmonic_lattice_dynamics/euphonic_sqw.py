
from euphonic import ForceConstants, ureg
from euphonic.util import mp_grid

import matplotlib.pyplot as plt
import numpy as np


def get_Qpts(Q_start,Q_end,Q_steps):

    num_Q_path = len(Q_start)
    num_Q_steps = sum(Q_steps)
    Qpts = np.zeros((num_Q_steps,3))

    s = 0
    for qq in range(num_Q_path):
        Qpts[s:s+Q_steps[qq],0] = np.linspace(Q_start[qq][0],Q_end[qq][0],Q_steps[qq]+1)[:-1]
        Qpts[s:s+Q_steps[qq],1] = np.linspace(Q_start[qq][1],Q_end[qq][1],Q_steps[qq]+1)[:-1]
        Qpts[s:s+Q_steps[qq],2] = np.linspace(Q_start[qq][2],Q_end[qq][2],Q_steps[qq]+1)[:-1]
        s += Q_steps[qq]

    return Qpts

# -----------------------------


# Qpts path 
Qpts = get_Qpts(Q_start=[[0.0,0.0,0.0],
                         [5.0,0.0,0.0],
                         [5.0,5.0,0.0],
                         [5.0,5.0,5.0]],
                Q_end=[[5.0,0.0,0.0],
                       [5.0,5.0,0.0],
                       [5.0,5.0,5.0],
                       [0.0,0.0,0.0]],
                Q_steps=[251,251,251,251])
nQ = Qpts.shape[0]

# energy axis
e_bins = np.arange(-100,100,0.1)*ureg('meV')
e_fwhm = 1.5*ureg('meV')
Q_fwhm = 0.05*ureg('1/angstrom')

# temperature
T = 300*ureg('K')

fc = ForceConstants.from_phonopy(path='./')

# get DW factor
q_grid = mp_grid([25,25,25])
dw = fc.calculate_qpoint_phonon_modes(q_grid,asr='reciprocal')
dw = dw.calculate_debye_waller(temperature=T)

# get phonon eigenvectors on path
sqw = fc.calculate_qpoint_phonon_modes(Qpts,asr='reciprocal')

# get SQW colormap
sqw = sqw.calculate_structure_factor(dw=dw)
sqw = sqw.calculate_sqw_map(e_bins,temperature=T)
sqw = sqw.broaden(y_width=e_fwhm)

# get the data from the object
sqw = sqw.to_dict()
e = sqw['y_data'][:]
sf = sqw['z_data'][:,:]

# plot it
fig, ax = plt.subplots(figsize=(12,6))

extent = (0,1,e.min(),e.max())
im = ax.imshow(sf.T,aspect='auto',origin='lower',cmap='hot',
        vmin=0,vmax=5,extent=extent)
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

plt.savefig('euphonic_sqw.pdf',bbox_inches='tight')






