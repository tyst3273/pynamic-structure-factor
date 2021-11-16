import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.io import loadmat
from matplotlib.image import imread
import sys

sys.path.append('/home/ty/python_modules/pynamic-structure-factor/')
from mod_io import read_sqw, read_bragg, read_timeavg

# ---------------------------------------------------------------------

def unfold_cut(in_file,num_Q,Q_range,data_type='timeavg',e_range=None,out_file=None):

    # check that the requested type is allowed
    if data_type not in ['bragg','sqw','timeavg']:
        print('\n ** ERROR **\n data type \'{data_type}\' is unknown\n')

    # get the data
    if data_type == 'sqw':
        energy, Qpts, raw = read_sqw(in_file)
    elif data_type == 'timeavg':
        Qpts, raw = read_timeavg(in_file)
    else:
        Qpts, raw = read_bragg(in_file)

    # if sqw, integrate (well, average) over energy
    if data_type == 'sqw':

        # use the full energy range if not requested. note, this is equivalent to the timeavg 
        # intensity, so it is a waste of time to do it this way.
        if e_range == None: 
            e_range = [energy.min(),energy.max()]

        # get the inds that fall inside the energy range
        e_inds = np.argwhere(energy >= e_range[0]).flatten()
        e_inds = np.intersect1d(np.argwhere(energy <= e_range[1]).flatten(),e_inds)

        # slice the data
        raw = raw[e_inds,:]
        raw = np.mean(raw,axis=0)

    # Q arrays used to assign data to elements in the unfolded array
    h_arr = np.linspace(Q_range[0],Q_range[1],num_Q)
    k_arr = np.linspace(Q_range[2],Q_range[3],num_Q)
    dQ = h_arr[1]-h_arr[0]
#    print('\ndQ = ',dQ,'\n')

    folded = np.zeros((num_Q,num_Q))
    for Qh in range(num_Q):
#        print(Qh,'/',num_Q)
        for Qk in range(num_Q):
            h = np.round(h_arr[Qh],4) # precision of the Q points in file
            k = np.round(k_arr[Qk],4) 
            hi = np.argwhere(Qpts[:,0]==h).flatten()
            ki = np.argwhere(Qpts[:,1]==k).flatten()
            ind = np.intersect1d(hi,ki).flatten()
            if ind.size == 0:
                continue
            else:
                folded[Qh,Qk] = raw[ind]

    # unfold the intensity from a right triangle onto positive quadrant in x-y plane
    tmp = np.copy(folded)
    np.fill_diagonal(tmp,0)
    folded = folded+tmp.T 

    # now unfold onto other quadrants (x,-y), (-x,y), and (-x,-y)
    unfolded = np.zeros((num_Q*2-1,num_Q*2-1))
    unfolded[:num_Q,:num_Q] = np.fliplr(np.flipud(folded[:,:]))
    unfolded[:num_Q,num_Q:] = np.flipud(folded[:,1:])
    unfolded[num_Q:,:num_Q] = np.fliplr(folded[1:,:])
    unfolded[num_Q:,num_Q:] = folded[1:,1:]

    # save the data to a txt file if requested
    if out_file != None:
        np.savetxt(out_file,unfolded,fmt='%2.6f',header=f'({-Q_range[1]:2.2f}, {Q_range[1]:2.2f}) '\
                '=> ({-Q_range[3]:2.2f}, {Q_range[3]:2.2f})')

    return unfolded






# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':


    # plot opts
    interp = 'none'
    cmap = 'viridis'


    # -------------------------------------------------------------------------------------
    # SQW files


    if False:

        # whether to unfold or load from txt
        unfold_new = False

        # opts to unfold/integrate
        num_Q = 121
        temp = 300
        HK = 6
        Q_range = [0,HK,0,HK]
        E = 0.5 # integrate from +- {E} meV

        # opts for plotting
        vlim_dict = {'ins_all'  :  [0,100],
                     'ins_cage' :  [0,20],
                     'ins_ma'   :  [0,30],
                     'xray'     :  [0,1500]}
        extent = [-Q_range[1],Q_range[1],-Q_range[3],Q_range[3]]


        cuts = ['deut_300/deut_300_ins_all_HK6.00_L0.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_ins_all_HK6.00_L1.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_ins_all_HK6.00_L2.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_ins_cage_HK6.00_L0.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_ins_cage_HK6.00_L1.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_ins_cage_HK6.00_L2.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_ins_ma_HK6.00_L0.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_ins_ma_HK6.00_L1.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_ins_ma_HK6.00_L2.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_xray_HK6.00_L0.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_xray_HK6.00_L1.50_SQW_FINAL.hdf5',
                'deut_300/deut_300_xray_HK6.00_L2.50_SQW_FINAL.hdf5']

        plot_file = f'pdfs/all_deut_300K_Epm_{E}meV_HK_6.pdf'
        pdf_doc = PdfPages(plot_file)

        for cut in cuts:

            f_name = cut.split('/')[1]
            print(f_name)

            which = f_name.strip('deut_300_').split('_HK')[0]
            L = f_name.split('_L')[1].strip('_SQW_FINAL.hdf5')
        
            out_file = f'unfolded_csv/deut_T_300K_{which}_Epm_{E}meV_HK_{HK}_L_{L}.txt'

            if unfold_new:
                unfolded = unfold_cut(cut,num_Q,Q_range,out_file=out_file,data_type='sqw',e_range=[-E,E])
            else:
                unfolded = np.loadtxt(out_file)

            extent = [-Q_range[1],Q_range[1],-Q_range[3],Q_range[3]] # extent of unfolded Q
            vlims = vlim_dict[which] 

            # set up and plot
            fig,ax=plt.subplots(figsize=(7.5,6))
            im = ax.imshow(unfolded,aspect='auto',cmap=cmap,origin='lower',extent=extent,
                    interpolation=interp,vmin=vlims[0],vmax=vlims[1])
            fig.colorbar(im,ax=ax,extend='both')

            # configure plot
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(1.5)
            ax.minorticks_on()
            ax.tick_params(which='both',width=1,labelsize='large')
            ax.tick_params(which='major',length=5)
            ax.tick_params(which='minor',length=2)

            xlabel = r'$\xi$ (r.l.u.)'
            ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')

            ylabel = r'$\eta$ (r.l.u.)'
            ax.set_ylabel(ylabel,labelpad=-2.0,fontweight='normal',fontsize='large')

            ax.annotate(which,xy=(-0.05,1.03),xycoords="axes fraction",color='b',fontsize='xx-large')
            fig.suptitle(rf"$\bf{{Q}}$=($\xi$,$\eta$,{L}), T={temp}K, -{E}<E<{E} meV",
                    fontsize='x-large',y=0.93)

            # save the plot
            pdf_doc.savefig(fig,dpi=150,bbox_inches='tight')

        # write this pdf
        pdf_doc.close()

        plt.close()

    # -------------------------------------------------------------------------------------
    # TIMEAVG files

    if False:

        # whether to unfold or load from txt
        unfold_new = False

        # opts to unfold/integrate
        num_Q = 121
        temp = 300
        HK = 6
        Q_range = [0,HK,0,HK]
        E = 5 # integrate from +- {E} meV

        # opts for plotting
        vlim_dict = {'ins_all'  :  [0,7.5],
                     'ins_cage' :  [0,1.5],
                     'ins_ma'   :  [0,6],
                     'xray'     :  [0,50]}
        extent = [-Q_range[1],Q_range[1],-Q_range[3],Q_range[3]]


        cuts = ['deut_300/deut_300_ins_all_HK6.00_L0.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_ins_all_HK6.00_L1.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_ins_all_HK6.00_L2.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_ins_cage_HK6.00_L0.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_ins_cage_HK6.00_L1.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_ins_cage_HK6.00_L2.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_ins_ma_HK6.00_L0.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_ins_ma_HK6.00_L1.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_ins_ma_HK6.00_L2.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_xray_HK6.00_L0.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_xray_HK6.00_L1.50_TIMEAVG_FINAL.hdf5',
                'deut_300/deut_300_xray_HK6.00_L2.50_TIMEAVG_FINAL.hdf5']

        plot_file = 'pdfs/all_deut_300K_diffuse_HK_6.pdf'
        pdf_doc = PdfPages(plot_file)

        for cut in cuts:

            f_name = cut.split('/')[1]
            print(f_name)

            which = f_name.strip('deut_300_').split('_HK')[0]
            L = f_name.split('_L')[1].strip('_TIMEAVG_FINAL.hdf5')

            out_file = f'unfolded_csv/deut_T_300K_{which}_diffuse_HK_{HK}_L_{L}.txt'

            if unfold_new:
                unfolded = unfold_cut(cut,num_Q,Q_range,out_file=out_file,data_type='timeavg')
            else:
                unfolded = np.loadtxt(out_file)

            extent = [-Q_range[1],Q_range[1],-Q_range[3],Q_range[3]] # extent of unfolded Q
            vlims = vlim_dict[which]

            # set up and plot
            fig,ax=plt.subplots(figsize=(7.5,6))
            im = ax.imshow(unfolded,aspect='auto',cmap=cmap,origin='lower',extent=extent,
                    interpolation=interp,vmin=vlims[0],vmax=vlims[1])
            fig.colorbar(im,ax=ax,extend='both')

            # configure plot
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(1.5)
            ax.minorticks_on()
            ax.tick_params(which='both',width=1,labelsize='large')
            ax.tick_params(which='major',length=5)
            ax.tick_params(which='minor',length=2)

            xlabel = r'$\xi$ (r.l.u.)'
            ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')

            ylabel = r'$\eta$ (r.l.u.)'
            ax.set_ylabel(ylabel,labelpad=-2.0,fontweight='normal',fontsize='large')

            ax.annotate(which,xy=(-0.05,1.03),xycoords="axes fraction",color='b',fontsize='xx-large')
            fig.suptitle(rf"$\bf{{Q}}$=($\xi$,$\eta$,{L}), T={temp}K, diffuse",
                    fontsize='x-large',y=0.93)

            # save the plot
            pdf_doc.savefig(fig,dpi=150,bbox_inches='tight')

            plt.close()

        # write this pdf
        pdf_doc.close()

        # -------------------------------------------------------------------------------------


    # -------------------------------------------------------------------------------------
    # ins-xray subtracted files

    if False:

        # opts to unfold/integrate
        num_Q = 121
        temp = 300
        HK = 6
        Q_range = [0,HK,0,HK]

        # opts for plotting
        xray_scale = 50
        ins_scale = 7.5
        which = 'ins-xray'
        vlim_dict = {'ins-xray' :  [0,1]}
        extent = [-Q_range[1],Q_range[1],-Q_range[3],Q_range[3]]


        ins_cuts =  ['unfolded_csv/deut_T_300K_ins_all_Epm_1meV_HK_6_L_0.50.txt',
                'unfolded_csv/deut_T_300K_ins_all_1meV_HK_6_L_1.50.txt',
                'unfolded_csv/deut_T_300K_ins_all_1meV_HK_6_L_2.50.txt']

        xray_cuts = ['unfolded_csv/deut_T_300K_xray_1meV_HK_6_L_0.50.txt',
                'unfolded_csv/deut_T_300K_xray_1meV_HK_6_L_1.50.txt',
                'unfolded_csv/deut_T_300K_xray_1meV_HK_6_L_2.50.txt']

        plot_file = 'pdfs/all_deut_300K_subtracted_HK_6.pdf'
        pdf_doc = PdfPages(plot_file)

        for ins, xray in zip(ins_cuts,xray_cuts):

            f_name = ins.split('/')[1]
            print(f_name)

            L = f_name.split('_L')[1].strip('_.txt')

            ins_data = np.loadtxt(ins)
            xray_data = np.loadtxt(xray)
            subtr = (ins_data/ins_data.max()-xray_data/xray_data.max())*40

            extent = [-Q_range[1],Q_range[1],-Q_range[3],Q_range[3]] # extent of unfolded Q
            vlims = vlim_dict[which]

            # set up and plot
            fig,ax=plt.subplots(figsize=(7.5,6))
            im = ax.imshow(subtr,aspect='auto',cmap=cmap,origin='lower',extent=extent,
                    interpolation=interp,vmin=vlims[0],vmax=vlims[1])
            fig.colorbar(im,ax=ax,extend='both')

            # configure plot
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(1.5)
            ax.minorticks_on()
            ax.tick_params(which='both',width=1,labelsize='large')
            ax.tick_params(which='major',length=5)
            ax.tick_params(which='minor',length=2)

            xlabel = r'$\xi$ (r.l.u.)'
            ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')

            ylabel = r'$\eta$ (r.l.u.)'
            ax.set_ylabel(ylabel,labelpad=-2.0,fontweight='normal',fontsize='large')

            ax.annotate(which,xy=(-0.05,1.03),xycoords="axes fraction",color='b',fontsize='xx-large')
            fig.suptitle(rf"$\bf{{Q}}$=($\xi$,$\eta$,{L}), T={temp}K, subtr",
                    fontsize='x-large',y=0.93)

            # save the plot
            pdf_doc.savefig(fig,dpi=150,bbox_inches='tight')

            plt.close()

        # write this pdf
        pdf_doc.close()

        # ----------------------------------------------------------------------------------------

        
    # -------------------------------------------------------------------------------------
    # tot-(cage+ma) subtracted files

    if True:

        HK = 6
        Q_range = [0,HK,0,HK]
        temp = 300

        which = 'tot-(cage+ma)'
        vlim_dict = {'tot-(cage+ma)' :  [0,1]}
        extent = [-Q_range[1],Q_range[1],-Q_range[3],Q_range[3]]

#        tot_cuts =  ['unfolded_csv/deut_T_300K_ins_all_diffuse_HK_6_L_0.50.txt',
#                'unfolded_csv/deut_T_300K_ins_all_diffuse_HK_6_L_1.50.txt',
#                'unfolded_csv/deut_T_300K_ins_all_diffuse_HK_6_L_2.50.txt']
#        ma_cuts = ['unfolded_csv/deut_T_300K_ins_ma_diffuse_HK_6_L_0.50.txt',
#                'unfolded_csv/deut_T_300K_ins_ma_diffuse_HK_6_L_1.50.txt',
#                'unfolded_csv/deut_T_300K_ins_ma_diffuse_HK_6_L_2.50.txt']
#        cage_cuts = ['unfolded_csv/deut_T_300K_ins_cage_diffuse_HK_6_L_0.50.txt',
#                'unfolded_csv/deut_T_300K_ins_cage_diffuse_HK_6_L_1.50.txt',
#                'unfolded_csv/deut_T_300K_ins_cage_diffuse_HK_6_L_2.50.txt']
#        plot_file = 'pdfs/all_deut_300K_diff_HK_6.pdf'
#        pdf_doc = PdfPages(plot_file)

        tot_cuts =  ['unfolded_csv/deut_T_300K_ins_all_Epm_1meV_HK_6_L_0.50.txt',
                'unfolded_csv/deut_T_300K_ins_all_Epm_1meV_HK_6_L_1.50.txt',
                'unfolded_csv/deut_T_300K_ins_all_Epm_1meV_HK_6_L_2.50.txt']
        ma_cuts = ['unfolded_csv/deut_T_300K_ins_ma_Epm_1meV_HK_6_L_0.50.txt',
                'unfolded_csv/deut_T_300K_ins_ma_Epm_1meV_HK_6_L_1.50.txt',
                'unfolded_csv/deut_T_300K_ins_ma_Epm_1meV_HK_6_L_2.50.txt']
        cage_cuts = ['unfolded_csv/deut_T_300K_ins_cage_Epm_1meV_HK_6_L_0.50.txt',
                'unfolded_csv/deut_T_300K_ins_cage_Epm_1meV_HK_6_L_1.50.txt',
                'unfolded_csv/deut_T_300K_ins_cage_Epm_1meV_HK_6_L_2.50.txt']
        plot_file = 'pdfs/all_deut_300K_diff_Epm_1meV_HK_6.pdf'
        pdf_doc = PdfPages(plot_file)

        for tot, ma, cage in zip(tot_cuts,ma_cuts,cage_cuts):

            f_name = tot.split('/')[1]
            print(f_name)

            L = f_name.split('_L')[1].strip('_.txt')

            tot_data = np.loadtxt(tot)
            ma_data = np.loadtxt(ma)
            cage_data = np.loadtxt(cage)
            diff = (tot_data-(ma_data+cage_data))

            np.savetxt(f'diff_{L}',diff)

            extent = [-Q_range[1],Q_range[1],-Q_range[3],Q_range[3]] # extent of unfolded Q
            vlims = vlim_dict[which]

            # set up and plot
            fig,ax=plt.subplots(figsize=(7.5,6))
            im = ax.imshow(diff,aspect='auto',cmap=cmap,origin='lower',extent=extent,
                    interpolation=interp,vmin=vlims[0],vmax=vlims[1])
            fig.colorbar(im,ax=ax,extend='both')

            # configure plot
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(1.5)
            ax.minorticks_on()
            ax.tick_params(which='both',width=1,labelsize='large')
            ax.tick_params(which='major',length=5)
            ax.tick_params(which='minor',length=2)

            xlabel = r'$\xi$ (r.l.u.)'
            ax.set_xlabel(xlabel,labelpad=4.0,fontweight='normal',fontsize='large')

            ylabel = r'$\eta$ (r.l.u.)'
            ax.set_ylabel(ylabel,labelpad=-2.0,fontweight='normal',fontsize='large')
        
            note = r'$\rho_{ma}$($\bf{Q}$,$\omega$)$^*\rho_{cage}$($\bf{Q}$,$\omega$)+'\
                    r'$\rho_{ma}$($\bf{Q}$,$\omega$)$\rho_{cage}$($\bf{Q}$,$\omega$)$^*$'
#            ax.annotate(note,xy=(0.04,0.95),xycoords="axes fraction",color='k',fontsize='x-large')
            fig.suptitle(rf"$\bf{{Q}}$=($\xi$,$\eta$,{L}), T={temp}K"+'\n'+note,
                    fontsize='x-large',y=0.98)

            # save the plot
            pdf_doc.savefig(fig,dpi=150,bbox_inches='tight')

            plt.close()

        # write this pdf
        pdf_doc.close()

        # ----------------------------------------------------------------------------------------

