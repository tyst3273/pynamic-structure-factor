import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt


class plot_from_txt:

    def __init__(self,sqw,Qmin,Qmax,fig_name=False,vmin=False,
            vmax=False,cmap='jet',interp='gaussian'):

        max_freq = sqw[:,0].max()
        min_freq = sqw[:,0].min()
        num_freq = sqw.shape[0]
        sqw = sqw[:,1:]

        fig,ax = plt.subplots()
        fig.set_size_inches(4,4,forward=True)

        if vmin == False or vmax == False:
            im = ax.imshow(np.log(sqw),aspect='auto',cmap=cmap,
                    origin='lower',extent=[0,1,min_freq,max_freq],interpolation=interp)
        else:
            im = ax.imshow(np.log(sqw),aspect='auto',cmap=cmap,
                    origin='lower',extent=[0,1,min_freq,max_freq],interpolation=interp,
                    vmin=vmin,vmax=vmax)

        fig.colorbar(im,ax=ax,extend='both')

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        ax.minorticks_on()
        ax.tick_params(which='both',width=1,labelsize='large')
        ax.tick_params(which='major',length=5)
        ax.tick_params(which='minor',length=2)

#        ax.axis([0,1,-1,10])

        ax.set_xticklabels([])
        xlabel = r'$\bf{{Q}}_i$=({:.2f},{:.2f},{:.2f})$\rightarrow$ '.format(Qmin[0],Qmin[1],
                Qmin[2])+r'$\bf{{Q}}_f$=({:.2f},{:.2f},{:.2f})'.format(Qmax[0],Qmax[1],Qmax[2])
        ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')
        ax.set_ylabel('Energy (meV)',labelpad=2.0,fontweight='normal',fontsize='large')

        plt.title(r'log(S($\bf{Q}$,$\omega$))',fontsize='x-large')

        if fig_name != False:
            plt.savefig(fig_name+'.pdf',format='pdf',dpi=300,bbox_inches='tight')
        plt.show()



class plot_single_from_txt:

    def __init__(self,sqw,Qpoint,sigma=False,fig_name=False,log_scale=False):

        max_freq = sqw[:,0].max()
        min_freq = sqw[:,0].min()
        num_freq = sqw.shape[0]
        energy = sqw[:,0]
        sqw = sqw[:,1:]
        if sigma != False:
            sqw = gaussian_filter(sqw,sigma)

        fig,ax = plt.subplots()
        fig.set_size_inches(6,4,forward=True)

        if log_scale == False:
            ax.plot(energy,sqw[:,0],color='k',lw=2,ls='-')
        else:
            ax.semilogy(energy,sqw[:,0],color='k',lw=2,ls='-')


        #####
        ax.plot([5.8,5.801],[0,30],color='r',ls='--')
        ax.plot([8,8.001],[0,30],color='b',ls='--')

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        ax.minorticks_on()
        ax.tick_params(which='both',width=1,labelsize='large')
        ax.tick_params(which='major',length=5)
        ax.tick_params(which='minor',length=2)
#        ax.set_yticklabels([])
#        ax.axis([1,25,4,20])
        
        ax.annotate(r'$\bf{{Q}}_i$=({:.2f},{:.2f},{:.2f})'.format(Qpoint[0],Qpoint[1],Qpoint[2]),
                xy=(0.5,0.85),xycoords="axes fraction",color='k',fontsize='x-large')

        ax.set_ylabel('Intensity (arb. units)',labelpad=4.0,
                                    fontweight='normal',fontsize='x-large')
        ax.set_xlabel('Energy (meV)',labelpad=2.0,fontweight='normal',fontsize='x-large')

        if fig_name != False:
            plt.savefig(fig_name+'.pdf',format='pdf',dpi=150,bbox_inches='tight')

        plt.show()

