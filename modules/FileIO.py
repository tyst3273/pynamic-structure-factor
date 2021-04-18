import numpy as np
import os

# ========================================================
# ------------ write result of calculation ---------------
# ========================================================

def save_sqw(params,Q_points,sqw,fmt='%6.10f',f_name='sqw.dat'):

    if not os.path.exists(params.output_dir):
        os.mkdir(params.output_dir)

    energy = params.meV

    header = ' 1st 3 rows are Qpts (the E col is (h,k,l)==(0,0,0)). the 1st col is E in meV'

    write_dat = np.zeros((energy.shape[0]+3,Q_points.shape[0]+1))
    write_dat[3:,0] = energy
    write_dat[:3,1:] = Q_points.T
    write_dat[3:,1:] = sqw

    f_name = os.path.join(params.output_dir,f_name)
    np.savetxt(f_name,write_dat,fmt=fmt,header=header)



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

    logo = r"""
          /##########\
         /###/  /#####\  /##/####|  /##/#####|/###||###\\############//########\
        /###/  /###/\##\/##/#####| /##/##/|##|##########|  |##|      /##/    \##\
       /##########/  \####/##/|##|/##/#######|##||##||##|  |##|     /##/
      /###/          /###/##/ |#####/##/  |##|##|    |##|  |##|    |##|
     /###/          /###/##/  |####/##/   |##|##|    |##|#######|  |###############\
    /###/#############################|   |##|##################|  |##|######\|#####\
   /###/##/       /###/    /##/ /#####|   |##|##|      /##/  |##|  |##|##| |##|##|
  /###/########/ /###/    /#######/|##|   |##|##|     /##/   |##|  |##|#######|#######\
 |###|     /##/ /###/    /##/ \##\ |##|   |##|##|    /##/    |##|  |##|##| \##\##|
 |###|#######/ /###/    /##/   \##\\########/ \####//##/      \######/ ##|  \###########|
 |#####################|####################\ /#######/      
 |###|       /###/  |##|##|     |##||##|  |##|##| |##|   |@@\    /@@| /@@@@|     /@@@@@\
 |#####################|##|     |##||##|  |##|#######|    \@@\  /@@/ /@/|@@|    |@@| |@@|
 |###|     /###/    |##|##|     |##||##|  |##|##| \##\     \@@\/@@/     |@@|    |@@| |@@|
 |###|    /###/     |##|\######||##| \######/|##|  \##\     \@@@@/   |@@@@@@@|(@)\@@@@@/

"""
    herald = """
 Pynamic Structure Factor, version 1.0

 Now with MPI support!

 Author: Ty Sterling
         Department of Physics
         University of Colorado Boulder

 Email: ty.sterling@colorado.edu
"""
    license = """
 This is free software licensed under the Gnu General Public License (GPL) version 3.
 You should have a copy of the license agreement in the top directory. 

 This software is presented in the hope that you will find it useful, but with no 
 warranty or gaurantee that the results will be accurate or valid. If you do happen 
 to notice bugs or errors, please notify me and I will try to address them (but I can 
 make no promises). Thanks and have fun :)

 """

    banner = """
 ########################################################################################
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ########################################################################################
"""
    print(banner)
    print(logo)
    print(herald)
    print(license)
    print(f' ** running the S(Q,w) calculation using {n_ranks} processes **\n')
    print(banner,flush=True)

    





