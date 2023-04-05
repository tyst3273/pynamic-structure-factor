
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"


from psf.m_io import c_reader


fig_name = 'si_xrays_phonons.pdf'

reader = c_reader('./../out/pristine_xrays_STRUFACS.hdf5')
reader.read_sqw()
Q_rlu = reader.Q_rlu
energy = reader.energy
pn_sqw = reader.sqw

reader = c_reader('./../out/substitutions_xrays_STRUFACS.hdf5')
reader.read_sqw()
Q_rlu = reader.Q_rlu
sn_sqw = reader.sqw


num_verts = 3
num_e = energy.size//2

cmap = 'viridis'
interp = 'none'
vmin = 0
vmax = 5e-5
c = (1,1,1)

fig, ax = plt.subplots(2,1,figsize=(7,7),gridspec_kw={'wspace':0.1,'hspace':0.1})


extent = [0,num_verts,0,energy.max()]


ax[0].imshow(pn_sqw[:,:num_e].T,origin='lower',aspect='auto',interpolation=interp,cmap=cmap,
                 vmin=vmin,vmax=vmax,extent=extent)
im = ax[1].imshow(sn_sqw[:,:num_e].T,origin='lower',aspect='auto',interpolation=interp,cmap=cmap,
                 vmin=vmin,vmax=vmax,extent=extent)

fig.colorbar(im,ax=[ax[0],ax[1]],location='right',extend='both',aspect=40,pad=0.03)


for ii in range(2):

    for axis in ['top','bottom','left','right']:
        ax[ii].spines[axis].set_linewidth(1.5)
    ax[ii].minorticks_on()
    ax[ii].tick_params(which='both',width=1,labelsize='medium')
    ax[ii].tick_params(which='major',length=5)
    ax[ii].tick_params(which='minor',length=2)
    ax[ii].set_ylim(0,70)
    ax[ii].set_xlim(0,3)
    ax[ii].set_xticks([0,1,2,3])
    ax[ii].set_rasterized(True)
    # ax[ii].set_ylim(1e-5,0.5)
    # ax[ii].set_yscale('log')


ax[0].plot([2,2],[0,energy.max()],lw=2,ls='-',c='k')
ax[1].plot([2,2],[0,energy.max()],lw=2,ls='-',c='k')

ax[1].set_xticklabels(['0,0,0','5,0,0','5,5,0; 2.5,0,0','2.5,5,0'])
ax[1].set_xlabel(r'$\bm{Q}$ (rlu)',labelpad=4,fontsize='large')

ax[0].set_xticklabels([])

ax[0].set_ylabel('Energy (meV)',labelpad=4,fontsize='large')
ax[1].set_ylabel('Energy (meV)',labelpad=4,fontsize='large')


ax[0].annotate(r'pristine',xy=(0.05,0.05),xycoords='axes fraction',fontsize='large',color=c)
ax[1].annotate('10\% Ge\nsubsitutions',xy=(0.025,0.05),
                           xycoords='axes fraction',fontsize='large',color=c)


plt.savefig(fig_name,dpi=150,bbox_inches='tight')


plt.show()

