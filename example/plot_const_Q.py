import sys
sys.path.append('/home/ty/research/repos/pynamic-structure-factor/modules')
import matplotlib.pyplot as plt
from mod_io import read_sqw

#axis_lims = [0,60,1e-6,1e2]


f_name = 'sqw/bands_SQW_FINAL.hdf5'
which_Q = 41

energy, Qpts, sqw = read_sqw(f_name)
num_e = energy.shape[0]//2
num_Q = Qpts.shape[0]
energy = energy[:num_e]
e_max = energy[-1]
sqw = sqw[:num_e,:]

Q = Qpts[which_Q,:] 

# create plots 
fig,ax = plt.subplots(figsize=(6,4))
ax.plot(energy,sqw[:,which_Q],marker='o',ms=2,mfc='b',mec='k',color='b',lw=1,ls=':')

# format plots 
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
#ax.axis(axis_lims)

xlabel = r'Energy (meV)'
ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
ylabel = 'intensity (arb. units.)'
ax.set_ylabel(ylabel,labelpad=2.0,fontweight='normal',fontsize='large')

#ax.annotate(r"all modes",xy=(0.03,0.925),
#        xycoords="axes fraction",color='w',fontsize='large')
# create plots 

fig.suptitle(rf"$\bf{{Q}}$=({Q[0]:2.2f},{Q[1]:2.2f},{Q[2]:2.2f})",fontsize='x-large',y=0.98)
plt.savefig('const_Q.pdf',format='pdf',dpi=150,bbox_inches='tight')

plt.show()



