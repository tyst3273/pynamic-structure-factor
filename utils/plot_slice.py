import sys
sys.path.append('./modules')
import matplotlib.pyplot as plt
from mod_io import read_sqw

vlims = [None,None]
vlims = [0,0.025]

axis_lims = [0,4,0,80]

interp = 'none'
cmap = 'CMRmap'


# load data 
if len(sys.argv) != 2:
    print('give file name to plot')
    exit()
else:
    f_name = sys.argv[1]


energy, Qpts, sqw = read_sqw(f_name)

num_e = energy.shape[0]//2
num_Q = Qpts.shape[0]

energy = energy[:num_e]
e_max = energy[-1]
sqw = sqw[:num_e,:]


# create plots 
fig,ax = plt.subplots(figsize=(4,4))
im = ax.imshow(sqw,aspect='auto',cmap=cmap,origin='lower',extent=[0,4,0,e_max],
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

xlabel = r'($\xi$,0,0) (r.l.u.)'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = r'Energy (meV)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')

#ax.annotate(r"all modes",xy=(0.03,0.925),
#        xycoords="axes fraction",color='w',fontsize='large')
#fig.suptitle(r"$\bf{Q}$=(0+$\xi$,0+$\xi$,0)",fontsize='x-large',y=0.98)
#plt.savefig('H0K0L0_H2K2L0_10K.pdf',format='pdf',dpi=150,bbox_inches='tight')

plt.show()



