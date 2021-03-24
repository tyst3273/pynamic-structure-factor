import numpy as np
import sys
sys.path.append('modules')
import Plot


sqw = np.loadtxt('sqw.dat')
with open('sqw.dat','r') as fid:
    tmp = fid.readline().strip().strip('#').split('...')
Qmin = tmp[0].strip().split()
Qmax = tmp[1].strip().split()

plot = Plot.plot_from_txt(sqw,Qmin,Qmax,fig_name='test')


