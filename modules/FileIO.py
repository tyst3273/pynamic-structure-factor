import numpy as np
import os

# ========================================================
# ------------ write result of calculation ---------------
# ========================================================

def save_sqw(params,sqw,fmt='%6.10f',f_name='sqw.dat'):

    if not os.path.exists(params.output_dir):
        os.mkdir(params.output_dir)

    header = 'meV '
    for Q in range(params.Qsteps):
        header = header+'{:2.3f} {:2.3f} {:2.3f} '.format(
            params.reduced_Q[Q,0],params.reduced_Q[Q,1],params.reduced_Q[Q,2])
    
    f_name = os.path.join(params.output_dir,f_name)
    np.savetxt(f_name,np.append(params.meV.reshape(params.num_freq,1),sqw,axis=1),fmt=fmt,header=header)



# ========================================================
# ------------ uniform way to print output ---------------
# ========================================================

def print_stdout(message,msg_type=None):

    """
    i actually print error messages with this too. kinda unforuntate naming, but too lazy
    to go fix it.
    """

    if msg_type != None:
        print(f'\n ** {msg_type} **')
    print(f' {message}',flush=True)



# ========================================================
# -------------- print herald on startup -----------------
# ========================================================

def print_herald(n_ranks):

    herald = """

 Pynamic Structure Factor, version 1.0

 Now with MPI support!

 Author: Ty Sterling
         Department of Physics
         University of Colorado Boulder

 Email: ty.sterling@colorado.edu
"""
    license = """

 This is free software licensed under the Gnu General Public License (GPL) v3. 
 You should have a copy of the license agreement in the top directory. 

 This software is presented in the hope that you will find it useful, but with 
 no warranty or gaurantee that the results will be accurate or valid. If you do
 happen to notice bugs or errors, please notify me and I will try to address 
 them (but I can make no promises). Thanks and have fun :)

 """

    banner = """

 #############################################################################
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #############################################################################
"""

    print(herald,flush=True)
    print(license,flush=True)
    print(f' running the S(Q,w) calculation using {n_ranks} processes',flush=True)
    print(banner,flush=True)





