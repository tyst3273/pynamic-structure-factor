from timeit import default_timer as timer
import os
import sys

path_to_modules = '/home/ty/custom_modules/pynamic-structure-factor/'
sys.path.append(path_to_modules)

import Parameters
import Calculator
import Plot
import FileIO




Qinds = [0]
temp = '10'
steps = [10] 
out_dir = f'sqe_cQ_{temp}K/'


# ===========    MD input parameters from calc   ==========

pos_dir = f'/mnt/5f27543f-dac3-492c-a181-1d9503ab021e/posvel/{temp}K/'
md_args = {'traj_file':pos_dir+f'pos_5_{temp}K.hdf5',
           'unwrap_pos':True,
           'recalculate_cell_lengths':True,
           'log_file':f'log.cQ_{temp}K_1',
           'dt':0.25e-15,
           'stride':200,
           'total_steps':4000000,
           'num_atoms':48000,
           'num_types':6,
           'supercell':[20,20,10],
           'b':{'1':6.646*1e-5,      # C  
                '2':9.360*1e-5,      # N
                '3':6.671*1e-5,      # H2, H = -3.741 
                '4':6.671*1e-5,      # H2, H = -3.741 
                '5':9.405*1e-5,      # Pb
                '6':5.280*1e-5},     # I
           'lattice_vectors':[[6.3,0,0], 
                              [0,6.3,0],
                              [0,0,6.3]]}



# ======  check if output dir exists, make it if not  ======= 

if not os.path.exists(out_dir):
    os.mkdir(out_dir)



# ================    loop over Q points    =====================

for Q in Qinds:

    sqw_args = {'Qpoints_file':f'Qi_{Q}', # if != False, overwrites Qmin, Qmax, Qsteps
                'blocks':200,
                'debug':True}


    # ===================    all modes  =====================

    md_args['b'] = {'1':6.646e-5,      # C  
                    '2':9.360e-5,      # N
                    '3':6.671e-5,      # H2 
                    '4':6.671e-5,      # H2 
                    '5':9.405e-5,      # Pb
                    '6':5.280e-5}     # I

    for step in steps:

        md_args['traj_file'] = pos_dir+f'pos_{step}_{temp}K.hdf5'
        params = Parameters.params(**md_args,**sqw_args) # initialize Qpoints, freq, etc.
        calc = Calculator.calc(params) # run the calculator. It is called within __init__() 
        
        calc.run(params)
        params.clean_up()

        f_name = (f"MAPI_{temp}_cQ_all_Q{Q}_{step}.dat")
        io = FileIO.io(params,calc,f_name=out_dir+f_name) # save SQE array as csv file


    # ===================   cage modes  =====================
 
    md_args['b'] = {'1':0,               # C
                    '2':0,               # N
                    '3':0,               # H2
                    '4':0,               # H2
                    '5':9.405e-5,        # Pb
                    '6':5.280e-5}        # I

    for step in steps:

        md_args['traj_file'] = pos_dir+f'pos_{step}_{temp}K.hdf5'
        params = Parameters.params(**md_args,**sqw_args) # initialize Qpoints, freq, etc.
        calc = Calculator.calc(params) # run the calculator. It is called within __init__()

        calc.run(params)
        params.clean_up()

        f_name = (f"MAPI_{temp}_cQ_cage_Q{Q}_{step}.dat")
        io = FileIO.io(params,calc,f_name=out_dir+f_name) # save SQE array as csv file


    # ===================   MA modes  =====================

    md_args['b'] = {'1':6.646e-5,      # C  
                    '2':9.360e-5,      # N
                    '3':6.671e-5,      # H2 
                    '4':6.671e-5,      # H2
                    '5':0,             # Pb
                    '6':0}             # I

    for step in steps:

        md_args['traj_file'] = pos_dir+f'pos_{step}_{temp}K.hdf5'
        params = Parameters.params(**md_args,**sqw_args) # initialize Qpoints, freq, etc.
        calc = Calculator.calc(params) # run the calculator. It is called within __init__()

        calc.run(params)
        params.clean_up()

        f_name = (f"MAPI_{temp}_cQ_ma_Q{Q}_{step}.dat")
        io = FileIO.io(params,calc,f_name=out_dir+f_name) # save SQE array as csv file






