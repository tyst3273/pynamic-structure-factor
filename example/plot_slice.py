import sys
sys.path.append('/home/ty/research/repos/pynamic-structure-factor/modules')
import matplotlib.pyplot as plt
from mod_io import read_sqw

#vlims = [None,None]
vlims = [0,0.00005]

interp = 'none'
cmap = 'CMRmap'

f_name = 'sqw/bands_SQW_FINAL.hdf5'
energy, Qpts, sqw = read_sqw(f_name)

num_e = energy.shape[0]//2
num_Q = Qpts.shape[0]

energy = energy[:num_e]
e_max = energy[-1]
sqw = sqw[:num_e,:]

axis_lims = [0,num_Q,0,e_max]

# create plots 
fig,ax = plt.subplots(figsize=(4,4))
im = ax.imshow(sqw,aspect='auto',cmap=cmap,origin='lower',extent=[0,num_Q,0,e_max],
        interpolation=interp,vmin=vlims[0],vmax=vlims[1])
fig.colorbar(im,ax=ax,extend='both')


# format plots 
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
ax.axis(axis_lims)

#xlabel = r'($\xi$,0,0) (r.l.u.)'
xlabel = r'$\bf{Q}$ index'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = r'Energy (meV)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')

plt.savefig('dispersions.pdf',format='pdf',dpi=150,bbox_inches='tight')

plt.show()



