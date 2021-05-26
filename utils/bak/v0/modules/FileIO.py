import numpy as np
import os

# ============ write result of calculation ===============

def save_sqe(params,sqw,fmt='%6.10f',f_name='sqw.dat'):

    if not os.path.exists(params.output_dir):
        os.mkdir(params.output_dir)

    header = 'meV '
    for Q in range(params.Qsteps):
        header = header+'{:2.3f} {:2.3f} {:2.3f} '.format(
            params.reduced_Q[Q,0],params.reduced_Q[Q,1],params.reduced_Q[Q,2])
    
    f_name = os.path.join(params.output_dir,f_name)
    np.savetxt(f_name,np.append(params.meV.reshape(params.num_freq,1),sqw,axis=1),fmt=fmt,header=header)



# ========== reading/writing to log file ===============

def print_log(log_handle,message,msg_type=None):

    if msg_type != None:
        log_handle.write(f'\n\n ** {msg_type} **\n')
    else:
        log_handle.write('\n')
    log_handle.write(f' {message}')
    log_handle.flush()


def print_herald(log_handle):

    herald = """

 Pynamic Structure Factor, version 1.0                 

 Author: Ty Sterling
         Department of Physics
         University of Colorado Boulder

 Email: ty.sterling@colorado.edu


 This is free software licensed under the Gnu General Public License (GPL) v3. 
 You should have a copy of the license agreement in the top directory. 

 This software is presented in the hope that you will find it useful, but with 
 no warranty or gaurantee that the results will be accurate or valid. If you do
 happen to notice bugs or errors, please notify me and I will try to address 
 them (but I can make no promises). Thanks and have fun :)

 #############################################################################
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #############################################################################

"""
    log_handle.write(herald)
    log_handle.flush()





