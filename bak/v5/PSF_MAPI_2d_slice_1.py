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
           'log_file':'log.2d',
           'dt':0.25e-15,
           'stride':200,
           'total_steps':4000000,
           'num_atoms':48000,
           'num_types':6,
           'supercell':[20,20,10],
           'b':{'1':0, #6.646*1e-5,  # C  # these are fm, converted to Angstrom later 
                '2':0, #9.360*1e-5,  # N
                '3':0, #6.671*1e-5,  # deuterated -3.741, # H
                '4':0, #6.671*1e-5,  # deuterated -3.741, # H
                '5':9.405*1e-5,  # Pb
                '6':5.280*1e-5}, # I
           'lattice_vectors':[[6.373,0,0], # using hdf5, these are read from the file
                              [0,6.373,0],
                              [0,0,6.373]]}

sqw_args = {'Qpoints_file':False,
            'Qmin':[0,0,1], 
            'Qmax':[0,2,1],
            'Qsteps':41,
            'blocks':20,
            'debug':False}


for step in [8,10]: #[5,6,7,8,9,10]:

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
    f_name = (f"sqe_MAPI_H{sqw_args['Qmin'][0]:2.2f}K{sqw_args['Qmin'][1]:2.2f}"
              f"L{sqw_args['Qmin'][2]:2.2f}_H{sqw_args['Qmax'][0]:2.2f}"
              f"K{sqw_args['Qmax'][1]:2.2f}L{sqw_args['Qmax'][2]:2.2f}_{step}_cage.dat")

    io = FileIO.io(params,calc,f_name=f'sqe/{temp}K/2d/'+f_name) # save SQE array as csv file





