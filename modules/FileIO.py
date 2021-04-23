import numpy as np
import os
import h5py

# ========================================================
# ------------ write result of calculation ---------------
# ========================================================

def save_sqw(params,Q_points,sqw,fmt='%6.10f',f_name='sqw.hdf5'):

    energy = params.meV
    num_e = params.num_freq
    num_Q = Q_points.shape[0]

    if not os.path.exists(params.output_dir):
        os.mkdir(params.output_dir)

    f_name = os.path.join(params.output_dir,f_name)
    with h5py.File(f_name,'w') as db:

        db_energy = db.create_dataset('energy_meV',[num_e])
        db_Qpts = db.create_dataset('Qpts_rlu',[num_Q,3])
        db_sqw = db.create_dataset('sqw',[num_e,num_Q])

        db_energy[:] = energy[:]
        db_Qpts[:,:] = Q_points[:,:]
        db_sqw[:,:] = sqw[:,:]


def read_sqw(f_name):
    
    with h5py.File(f_name,'r') as db:
        energy = db['energy_meV'][:]
        Qpts = db['Qpts_rlu'][:,:]
        sqw = db['sqw'][:,:]

    return energy, Qpts, sqw



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
#    print(banner)
#    print(logo)
    print(herald)
    print(license)
    print(f' ** running the S(Q,w) calculation using {n_ranks} processes **\n')
    print(banner,flush=True)

    





