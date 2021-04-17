import sys
import mpi4py

sys.path.append('modules')

import Parser
import Parameters 
import FileIO
import Calculator


print('\n starting the S(Q,w) calculation. Check the log file for progress.\n')


# parser input file. should be done on rank 1.
parser = Parser.parser()
parser.parse('input_params')


# initialize params object. should be done on each rank.
# should scatter array of Qpoints to each rank and run calc over those Q
params = Parameters.params(parser)

# initialize calculator
calc = Calculator.calc(params)
calc.run(params)


# close the files, etc
params.clean_up()



