
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"

from psf.m_io import c_reader


fig_name = 'rutile_total_neutrons.pdf'


reader = c_reader('./../out/pristine_neutrons_STRUFACS.hdf5')
reader.read_sqw()

H = reader.H
K = reader.K

energy = reader.energy
pn_sqw = reader.sqw[:,:,0,:]

pn_sqw = pn_sqw.sum(axis=2)


reader = c_reader('./../out/vacancies_neutrons_STRUFACS.hdf5')
reader.read_sqw()
vn_sqw = reader.sqw[:,:,0,:]

vn_sqw = vn_sqw.sum(axis=2)




num_verts = 3
num_e = energy.size//2

cmap = 'viridis'
interp = 'none'
vmin = 0
vmax = 0.3
c = (1,1,1)

fig, ax = plt.subplots(1,2,figsize=(8,3.25),gridspec_kw={'wspace':0.1,'hspace':0.1})


extent = [H.min(),H.max(),K.min(),K.max()]


ax[0].imshow(pn_sqw.T,origin='lower',aspect='auto',interpolation=interp,cmap=cmap,
                 vmin=vmin,vmax=vmax,extent=extent)
im = ax[1].imshow(vn_sqw.T,origin='lower',aspect='auto',interpolation=interp,cmap=cmap,
                 vmin=vmin,vmax=vmax,extent=extent)

fig.colorbar(im,ax=[ax[0],ax[1]],location='right',extend='both',aspect=40,pad=0.03)


for ii in range(2):

    for axis in ['top','bottom','left','right']:
        ax[ii].spines[axis].set_linewidth(1.5)
    ax[ii].minorticks_on()
    ax[ii].tick_params(which='both',width=1,labelsize='medium')
    ax[ii].tick_params(which='major',length=5)
    ax[ii].tick_params(which='minor',length=2)
    ax[ii].set_xticks([-3,-2,-1,0,1,2,3])
    ax[ii].set_yticks([-3,-2,-1,0,1,2,3])
    ax[ii].set_rasterized(True)


ax[1].set_yticklabels([])

ax[0].set_ylabel('H (rlu)',labelpad=4,fontsize='large')
ax[0].set_xlabel('K (rlu)',labelpad=4,fontsize='large')
ax[1].set_xlabel('K (rlu)',labelpad=4,fontsize='large')

ax[0].annotate(r'pristine',xy=(0.05,0.05),xycoords='axes fraction',fontsize='large',color=c)
ax[1].annotate('10\% O\nvacancies',xy=(0.02,0.02),
                           xycoords='axes fraction',fontsize='large',color=c)


plt.savefig(fig_name,dpi=150,bbox_inches='tight')


plt.show()

