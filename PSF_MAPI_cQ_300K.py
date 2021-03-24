from timeit import default_timer as timer
import os
import sys
sys.path.append('modules')
import Parameters
import Calculator
import Plot
import FileIO

# loops over files with lists of Q points, 
# computes S(Q,w) at each Q point in list,
# loops over 'steps', i.e. different MD 
# trajectory files, and saves output for 
# each Q-point list/step to .csv file in
# specified output directory.

# Q list indicies
Qinds = [0,1,2,3,4,5,6] #[7,8,9,10,11,12]
# temp to be calculated
temp = '300'
# MD files
steps = [8,10] 

# path to *.hdf5 files
pos_dir = f'/mnt/5f27543f-dac3-492c-a181-1d9503ab021e/posvel/{temp}K/'
# MD params needed by code
md_args = {'traj_file':pos_dir+f'pos_5_{temp}K.hdf5', # name of file
           'unwrap_pos':True,                         # un-impose minumum image if needed
           'recalculate_box_lengths':False,           # calculate box-lengths from MD calc
           'file_format':'hdf5',                      # not-needed anymore
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

# create output directory if it doesnt exist
if not os.path.exists(f'sqe_{temp}K'):
    os.mkdir(f'sqe_{temp}K')


# ===============================================================
# **************       loop over Q points      ******************
# ===============================================================

for Q in Qinds:

    sqw_args = {'Qpoints_file':f'Q_inputs/Qi_{Q}', # if != False, overwrites Qmin, Qmax, Qsteps
                'blocks':20,
                'debug':False}

    # =========================================================
    # **************      All_modes      **********************
    # =========================================================

    for step in steps:

        md_args['traj_file'] = pos_dir+f'pos_{step}_{temp}K.hdf5'
        params = Parameters.params(**md_args,**sqw_args) # initialize Qpoints, freq, etc.
        calc = Calculator.calc(params) # run the calculator. It is called within __init__() 
        f_name = (f"MAPI_{temp}_cQ_all_Q{Q}_{step}.dat")
        io = FileIO.io(params,calc,f_name=f'sqe_{temp}K/'+f_name) # save SQE array as csv file

    # =========================================================
    # **************      Cage_modes      *********************
    # =========================================================
    
    md_args['b'] = {'1':0,               # C
                    '2':0,               # N
                    '3':0,               # H2, H = -3.741 
                    '4':0,               # H2, H = -3.741
                    '5':9.405*1e-5,      # Pb
                    '6':5.280*1e-5}      # I

    for step in steps:

        md_args['traj_file'] = pos_dir+f'pos_{step}_{temp}K.hdf5'
        params = Parameters.params(**md_args,**sqw_args) # initialize Qpoints, freq, etc.
        calc = Calculator.calc(params) # run the calculator. It is called within __init__()
        f_name = (f"MAPI_{temp}_cQ_cage_Q{Q}_{step}.dat")
        io = FileIO.io(params,calc,f_name=f'sqe_{temp}K/'+f_name) # save SQE array as csv file

    # =========================================================
    # **************      MA_modes      ***********************
    # =========================================================

    md_args['b'] = {'1':6.646*1e-5,      # C  
                    '2':9.360*1e-5,      # N
                    '3':6.671*1e-5,      # H2, H = -3.741 
                    '4':6.671*1e-5,      # H2, H = -3.741 
                    '5':0,               # Pb
                    '6':0}               # I

    for step in steps:

        md_args['traj_file'] = pos_dir+f'pos_{step}_{temp}K.hdf5'
        params = Parameters.params(**md_args,**sqw_args) # initialize Qpoints, freq, etc.
        calc = Calculator.calc(params) # run the calculator. It is called within __init__()
        f_name = (f"MAPI_{temp}_cQ_MA_Q{Q}_{step}.dat")
        io = FileIO.io(params,calc,f_name=f'sqe_{temp}K/'+f_name) # save SQE array as csv file






