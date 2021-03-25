from timeit import default_timer as timer
import os
import sys

path_to_modules = '/home/ty/custom_modules/pynamic-structure-factor/'
sys.path.append(path_to_modules)

import Parameters
import Calculator
import Plot
import FileIO


temp = '300'
steps = [10]


# ===========    MD input parameters from calc   ==========

pos_dir = f'/mnt/5f27543f-dac3-492c-a181-1d9503ab021e/posvel/{temp}K/'
md_args = {'traj_file':pos_dir+f'pos_5_{temp}K.hdf5', # this gets over written for each step
           'unwrap_pos':True,
           'recalculate_cell_lengths':True,
           'log_file':f'log.2d',
           'dt':0.25e-15,
           'stride':200,
           'total_steps':4000000,
           'num_atoms':48000,
           'num_types':6,
           'supercell':[20,20,10],
           'b':{'1':6.646e-5,             # C  
                '2':9.360e-5,             # N
                '3':6.671e-5,             # H2
                '4':6.671e-5,             # H2 
                '5':9.405e-5,             # Pb
                '6':5.280e-5},            # I
           'lattice_vectors':[[6.3,0,0],
                              [0,6.3,0],
                              [0,0,6.3]]}


# =============   SQW calcultion options   ===============

sqw_args = {'Qmin':[0,0,0],
            'Qmax':[1,1,0],
            'Qsteps':21,
            'blocks':40,
            'debug':True}


# ======  check if output dir exists, make it if not  ======= 

if not os.path.exists(f'sqe_2d_{temp}K'):
    os.mkdir(f'sqe_2d_{temp}K')        


# ================    run the calculation   =================

for step in steps:

    md_args['traj_file'] = pos_dir+f'pos_{step}_{temp}K.hdf5'   # which file to use
    params = Parameters.params(**md_args,**sqw_args)            # initialize parameters
    calc = Calculator.calc(params)                              # initialize calculator

    calc.run(params)                                            # run calc
    params.clean_up()                                           # close files, probably not necessary

    f_name = (f"MAPI_{temp}_2d_all_{step}.dat")                 # output file name
    io = FileIO.io(params,calc,f_name=f'sqe_2d_{temp}K/'+f_name)   # save SQE array as csv file





