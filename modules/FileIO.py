import numpy as np

class io:

    def __init__(self,params,calc,fmt='%6.10f',f_name='sqw.dat'):

        """
        just save the SQW array as a *.csv file.
        """

        header = 'meV '
        for Q in range(params.Qsteps):
            header = header + '{:2.3f} {:2.3f} {:2.3f} '.format(params.reduced_Q[Q,0],
                    params.reduced_Q[Q,1],params.reduced_Q[Q,2])
        
        np.savetxt(f_name,np.append(params.meV.reshape(params.num_freq,1),
            calc.sqw,axis=1),fmt=fmt,header=header)
