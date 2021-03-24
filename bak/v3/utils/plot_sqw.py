import numpy as np
import sys
sys.path.append('modules')
import Plot


if len(sys.argv) != 1:
    f_name = sys.argv[1]
else:
    f_name = 'sqe.dat'

sqw = np.loadtxt(f_name)
with open(f_name,'r') as fid:
    tmp = fid.readline().strip().strip('#').split('...')
Qmin = tmp[0].strip().split()
Qmax = tmp[1].strip().split()

if sqw.shape[1] != 2:
    plot = Plot.plot_from_txt(sqw,Qmin,Qmax,fig_name='test',mask_tol=1e-6,mask_delta=1e-6)

else:
    plot = Plot.plot_single_from_txt(sqw,Qmin,fig_name='test_single',log_scale=False)


