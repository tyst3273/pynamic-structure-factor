import numpy as np
import matplotlib.pyplot as plt


class color_map:

    def __init__(self,params,calc,fig_name=False,mask_tol=1e-8):

        fig,ax = plt.subplots()
        fig.set_size_inches(6,6,forward=True)
        
        sqw = np.copy(calc.sqw)
        sqw = sqw*(sqw > 0).astype(int)
        sqw = sqw+(sqw <= mask_tol).astype(int)*mask_tol

        im = ax.imshow(np.log(sqw[:sqw.shape[0]//2,:]),aspect='auto',cmap='jet',
                origin='lower',extent=[0,1,0,params.max_freq],interpolation='gaussian')
        fig.colorbar(im,ax=ax,extend='both')

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        ax.minorticks_on()
        ax.tick_params(which='both',width=1,labelsize='large')
        ax.tick_params(which='major',length=5)
        ax.tick_params(which='minor',length=2)
        ax.set_xticklabels([])
        xlabel = r'$\bf{{Q}}_i$=({},{},{})$\rightarrow$ '.format(params.Qmin[0],
            params.Qmin[1],params.Qmin[2])+r'$\bf{{Q}}_f$=({},{},{})'.format(params.Qmax[0],
                    params.Qmax[1],params.Qmax[2])
        ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='x-large')
        ax.set_ylabel('Energy (meV)',labelpad=2.0,fontweight='normal',fontsize='x-large')

        if fig_name != False:
            plt.savefig(fig_name,format='png',dpi=150,bbox_inches='tight')

        plt.show()


class plot_from_txt:

    def __init__(self,sqw,Qmin,Qmax,fig_name=False):

        max_freq = sqw[:,0].max()/2
        sqw = sqw[:,1:]
        sqw = sqw+(sqw <= 1e-10).astype(int)*1e-10
        
        fig,ax = plt.subplots()
        fig.set_size_inches(8,6,forward=True)

        im = ax.imshow(np.log(sqw[:sqw.shape[0]//2,:]),aspect='auto',cmap='jet',
                origin='lower',extent=[0,1,0,max_freq],interpolation='gaussian')
        fig.colorbar(im,ax=ax,extend='both')

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        ax.minorticks_on()
        ax.tick_params(which='both',width=1,labelsize='large')
        ax.tick_params(which='major',length=5)
        ax.tick_params(which='minor',length=2)
        ax.set_xticklabels([])
        xlabel = r'$\bf{{Q}}_i$=({},{},{})$\rightarrow$ '.format(Qmin[0],Qmin[1],
                Qmin[2])+r'$\bf{{Q}}_f$=({},{},{})'.format(Qmax[0],Qmax[1],Qmax[2])
        ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='x-large')
        ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='x-large')
        ax.set_ylabel('Energy (meV)',labelpad=2.0,fontweight='normal',fontsize='x-large')

        if fig_name != False:
            plt.savefig(fig_name+'.png',format='png',dpi=150,bbox_inches='tight')

        plt.show()


