from timeit import default_timer as timer
import os
import sys
import numpy as np

path_to_modules = '/home/ty/custom_modules/pynamic-structure-factor/'
sys.path.append(path_to_modules)

import Parameters
import Calculator
import Plot
import FileIO



temp = '300'
step = 10
out_dir = f'sqe_2d_{temp}K_L_1.5'


pos_dir = f'/mnt/5f27543f-dac3-492c-a181-1d9503ab021e/posvel/{temp}K/'
md_args = {'traj_file':pos_dir+f'pos_{step}_{temp}K.hdf5',
           'dt':0.25e-15,
           'stride':200,
           'total_steps':4000000,
           'num_atoms':48000,
           'num_types':6,
           'supercell':[20,20,10]}

sqw_args = {'Qsteps':81,
            'num_blocks':20,
            'blocks':[0,2,4]}




if not os.path.exists(out_dir):
    os.mkdir(out_dir)        

for Qk in np.arange(0,4.05,0.05):

    Qk = np.round(Qk,2)
    print(f'\n ** Qi=(0,{Qk},1.5) -> Qf=(4,{Qk},1.5) **\n')


    sqw_args['Qmin'] = [0,Qk,1.5]
    sqw_args['Qmax'] = [4,Qk,1.5]


    # ================ all atoms ====================
    print(' ** All atoms ** ')

    md_args['b'] = {'1':6.646e-5,   # C
                    '2':9.360e-5,   # N
                    '3':6.671e-5,   # H2
                    '4':6.671e-5,   # H2
                    '5':9.405e-5,   # Pb
                    '6':5.280e-5}   # I

    params = Parameters.params(**md_args,**sqw_args)    
    calc = Calculator.calc(params) 

    calc.run(params)         
    params.clean_up()    

    f_name = f"MAPI_all_{temp}K_{step}_Qk_{Qk}.dat"
    io = FileIO.io(params,calc,
            f_name=os.path.join(out_dir,f_name)) 
   

    # ================ cage atoms ===================
    print(' ** Cage atoms ** ')

    md_args['b'] = {'1':0,          # C
                    '2':0,          # N
                    '3':0,          # H2
                    '4':0,          # H2
                    '5':9.405e-5,   # Pb
                    '6':5.280e-5}   # I

    params = Parameters.params(**md_args,**sqw_args)
    calc = Calculator.calc(params)

    calc.run(params)
    params.clean_up()

    f_name = f"MAPI_cage_{temp}K_{step}_Qk_{Qk}.dat"
    io = FileIO.io(params,calc,
            f_name=os.path.join(out_dir,f_name))


    # ================== ma atoms ===================
    print(' ** MA atoms ** ')

    md_args['b'] = {'1':6.646e-5,   # C
                    '2':9.360e-5,   # N
                    '3':6.671e-5,   # H2
                    '4':6.671e-5,   # H2
                    '5':0,          # Pb
                    '6':0}          # I

    params = Parameters.params(**md_args,**sqw_args)
    calc = Calculator.calc(params)

    calc.run(params)
    params.clean_up()

    f_name = f"MAPI_ma_{temp}K_{step}_Qk_{Qk}.dat"
    io = FileIO.io(params,calc,
            f_name=os.path.join(out_dir,f_name))








