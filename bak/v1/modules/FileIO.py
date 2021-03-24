import numpy as np

class io:

    def __init__(self,params,calc,fmt='%6.10f',f_name='sqw.dat'):
        
        np.savetxt(f_name,np.append(params.meV.reshape(params.num_freq,1),
            calc.sqw,axis=1),fmt=fmt,header='{} {} {} ... {} {} {}\n'.format(
                params.Qmin[0],params.Qmin[1],params.Qmin[2],params.Qmax[0],
                params.Qmax[1],params.Qmax[2])+'meV, S(Q_0,w), S(Q_1,w), S(Q_2,w), ...')
