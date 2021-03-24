from timeit import default_timer as timer
import os
import sys
sys.path.append('modules')
import Parameters
import Calculator
import Plot
import FileIO



##################################################################################
########################## INPUT PARAMETERS ######################################
##################################################################################

# these are the parameters during the simulation
# dt is in seconds, b and lattice vectors in angstrom
md_args = {'traj_file':'pos/pos_6_300K.hdf5',
           'file_format':'hdf5',
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

# these are parameters used to compute S(Q,w)
# Qmin, Qmax are in reciprocal lattice units
# Qstep is the number of steps along the BZ slice
# to do a single Qpoint, set Qsteps = 1 and Qmin to the Qpoint 
sqw_args = {'Qmin':[6,0,0],
            'Qmax':[6.5,0,0],
            'Qsteps':20,
            'blocks':20,
            'debug':True}



##################################################################################
############################# THE CALCULATION ####################################
##################################################################################

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




##################################################################################
############################## WRITE OUTPUT ######################################
##################################################################################

# create output directory if it doesnt exist
if not os.path.exists('sqe'):
    os.mkdir('sqe')

# output file name based on Qpoints 
f_name = (f"sqe_MAPI_H{sqw_args['Qmin'][0]:2.2f}K{sqw_args['Qmin'][1]:2.2f}"
          f"L{sqw_args['Qmin'][2]:2.2f}_H{sqw_args['Qmax'][0]:2.2f}"
          f"K{sqw_args['Qmax'][1]:2.2f}L{sqw_args['Qmax'][2]:2.2f}.dat")
io = FileIO.io(params,calc,f_name='sqe/'+f_name) # save SQE array as csv file





