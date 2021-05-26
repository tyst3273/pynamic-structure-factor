import sys
import numpy as np
from timeit import default_timer as timer

sys.path.append('/home/ty/custom_modules/pynamic-structure-factor/')

import Parser
import Parameters 
import FileIO
import Calculator
from ParalUtils import *


# get input file if given as arg
if len(sys.argv) != 1:
    input_file = sys.argv[1]
else:
    input_file = 'input_params'


# start a timer
start_time = timer()


# print the herald 
FileIO.print_herald(n_ranks=11)


# get input opts from file and initialize params
parser = Parser.parser()
parser.parse(input_file)
params = Parameters.params(parser)
params.rank_0_init()
total_reduced_Q = params.total_reduced_Q
params.rank_x_init(total_reduced_Q,rank=0)


# initiliaze and run the calculator
calc = Calculator.calc(params)
calc.run(params)


# save it the SQW result
f_name = params.outfile_prefix+f'_FINAL.hdf5'
FileIO.save_sqw(params,calc.sqw,f_name=f_name)


# calculate and print elapsed time
end_time = timer()
elapsed_time = (end_time-start_time)/60 #minutes
message = f'total elapsed time: {elapsed_time:2.3f} minutes'
FileIO.print_stdout(message,msg_type='TIMING')
 









