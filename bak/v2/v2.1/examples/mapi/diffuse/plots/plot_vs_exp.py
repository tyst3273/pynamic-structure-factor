
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"

from psf.m_io import c_reader

fig_name = 'mapi_vs_xrays.pdf'


reader = c_reader('./../out/protonated_xrays_STRUFACS.hdf5')
reader.read_energy_integrated_sqw(e_min=-2,e_max=2)

H = reader.H
K = reader.K
L = reader.L[0]

sq_integrated = reader.sq_integrated[:,:,0]


exp = np.loadtxt('xds_cut_L2.5')

cmap = 'viridis'
interp = 'none'
vmin = 0
vmax = 1.2
c = (1,1,1)

fig, ax = plt.subplots(1,2,figsize=(8,3.25),gridspec_kw={'wspace':0.1,'hspace':0.1})


extent = [H.min(),H.max(),K.min(),K.max()]


ax[0].imshow(sq_integrated.T,origin='lower',aspect='auto',interpolation=interp,cmap=cmap,
                 vmin=vmin,vmax=vmax,extent=extent)

vmax = 500
im = ax[1].imshow(exp.T,origin='lower',aspect='auto',interpolation=interp,cmap=cmap,
                 vmin=vmin,vmax=vmax,extent=extent)

fig.colorbar(im,ax=[ax[0],ax[1]],location='right',extend='both',aspect=40,pad=0.03)


for ii in range(2):

    for axis in ['top','bottom','left','right']:
        ax[ii].spines[axis].set_linewidth(1.5)
    ax[ii].minorticks_on()
    ax[ii].tick_params(which='both',width=1,labelsize='medium')
    ax[ii].tick_params(which='major',length=5)
    ax[ii].tick_params(which='minor',length=2)
    ax[ii].set_xticks([-4,-3,-2,-1,0,1,2,3,4])
    ax[ii].set_yticks([-4,-3,-2,-1,0,1,2,3,4])
    ax[ii].set_rasterized(True)


ax[1].set_yticklabels([])

ax[0].set_ylabel('K (rlu)',labelpad=4,fontsize='large')
ax[0].set_xlabel('H (rlu)',labelpad=4,fontsize='large')
ax[1].set_xlabel('H (rlu)',labelpad=4,fontsize='large')

ax[0].annotate(r'calculation',xy=(0.05,0.1),xycoords='axes fraction',fontsize='large',color=c)
ax[1].annotate('experiment',xy=(0.05,0.1),
                           xycoords='axes fraction',fontsize='large',color=c)

fig.suptitle(f'L={L:3.2f} (rlu)',fontsize='large')


plt.savefig(fig_name,dpi=150,bbox_inches='tight')

plt.show()

