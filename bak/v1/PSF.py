import sys
sys.path.append('modules')

import Parameters
import Calculator
import Plot
import FileIO


# these are the parameters during the simulation
# dt is in seconds, b and lattice vectors in angstrom
md_args = {'traj_file':'lammps/pos.dat',
           'dt':0.75e-15,
           'stride':32,
           'total_steps':2**21,
           'num_atoms':512,
           'num_types':2,
           'supercell':[16,2,2],
           'b':{'1':4.1491e-5},
           'lattice_vectors':[[5.431,0,0],
                              [0,5.431,0],
                              [0,0,5.431]]}

# these are parameters used to compute S(Q,w)
# Qmin, Qmax are in reciprocal lattice units
# Qstep is the number of steps along the BZ slice
sqw_args = {'Qmin':[2,0,0],
            'Qmax':[6,0,0],
            'Qsteps':128,
            'blocks':10,
            'debug':False}

params = Parameters.params(**md_args,**sqw_args)
calc = Calculator.calc(params)
io = FileIO.io(params,calc,f_name='sqw_1.dat')


