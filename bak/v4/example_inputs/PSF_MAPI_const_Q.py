from timeit import default_timer as timer
import os
import sys
sys.path.append('modules')
import Parameters
import Calculator
import Plot
import FileIO


temp = '300'
pos_dir = f'/mnt/5f27543f-dac3-492c-a181-1d9503ab021e/posvel/{temp}K/'

md_args = {'traj_file':pos_dir+f'pos_5_{temp}K.hdf5',
           'file_format':'hdf5',
           'log_file':'log.constQ',
           'dt':0.25e-15,
           'stride':200,
           'total_steps':4000000,
           'num_atoms':48000,
           'num_types':6,
           'supercell':[20,20,10],
           'b':{'1':6.646*1e-5,  # C  # these are fm, pos (Q) are in Angst (1/Angst) ... 
                '2':9.360*1e-5,  # N
                '3':6.671*1e-5,  # deuterated -3.741, # H
                '4':6.671*1e-5,  # deuterated -3.741, # H
                '5':9.405*1e-5,  # Pb
                '6':5.280*1e-5}, # I
           'lattice_vectors':[[6.373,0,0], # using hdf5, these are read from the file
                              [0,6.373,0],
                              [0,0,6.373]]}

sqw_args = {'Qpoints_file':'Qpoints', # if != False, overwrites Qmin, Qmax, Qsteps
            'blocks':10,
            'debug':False}


for step in [5,6,7,8,9,10]:

    md_args['traj_file'] = pos_dir+f'pos_{step}_{temp}K.hdf5'

    # start a timer
    start_time = timer()

    # do the calculation
    params = Parameters.params(**md_args,**sqw_args) # initialize Qpoints, freq, etc.
    calc = Calculator.calc(params) # run the calculator. It is called within __init__() 

    # print the elapsed time
    end_time = timer()
    elapsed_time = (end_time-start_time)/60 # minutes
    print(f'\n\tElapsed time:\t{elapsed_time:2.3f} minutes\n')
    params.log_handle.write(f'\n** Total elapsed time **\n {elapsed_time:2.3f} minutes **\n')
    params.log_handle.flush()

    # create output directory if it doesnt exist
    if not os.path.exists('sqe'):
        os.mkdir('sqe')

    # output file name based on Qpoints 
    f_name = (f"sqe_MAPI_constQ_{step}.dat")

    io = FileIO.io(params,calc,f_name=f'sqe/{temp}K/constQ/'+f_name) # save SQE array as csv file





