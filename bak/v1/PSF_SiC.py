import sys
sys.path.append('modules')

import Parameters
import Calculator
import Plot
import FileIO

# these are the parameters during the simulation
# dt is in seconds, b and lattice vectors in angstrom
md_args = {'traj_file':'lammps/SiC/2/pos.dat',
           'dt':0.5e-15,
           'stride':16,
           'total_steps':2**19,
           'num_atoms':2744,
           'num_types':2,
           'supercell':[7,7,7],
           'b':{'1':4.1491,'2':6.6460},
           'lattice_vectors':[[4.326,0,0],
                              [0,4.326,0],
                              [0,0,4.326]]}

# these are parameters used to compute S(Q,w)
# Qmin, Qmax are in reciprocal lattice units
# Qstep is the number of steps along the BZ slice
sqw_args = {'Qmin':[2,0,0],
            'Qmax':[6,0,0],
            'Qsteps':64,
            'blocks':20,
            'debug':True}

params = Parameters.params(**md_args,**sqw_args)
calc = Calculator.calc(params)
#plot = Plot.color_map(params,calc,mask_tol=1e-8)
io = FileIO.io(params,calc,f_name='sqe_7x7x7_1.dat')

sqw_args['Qmin'] = [2,2,0]
sqw_args['Qmax'] = [6,6,0]
params = Parameters.params(**md_args,**sqw_args)
calc = Calculator.calc(params)
io = FileIO.io(params,calc,f_name='sqe_7x7x7_2.dat')

sqw_args['Qmin'] = [2,2,2]
sqw_args['Qmax'] = [6,6,6]
params = Parameters.params(**md_args,**sqw_args)
calc = Calculator.calc(params)
io = FileIO.io(params,calc,f_name='sqe_7x7x7_3.dat')

sqw_args['Qmin'] = [2,0,1]
sqw_args['Qmax'] = [6,0,1]
params = Parameters.params(**md_args,**sqw_args)
calc = Calculator.calc(params)
io = FileIO.io(params,calc,f_name='sqe_7x7x7_4.dat')

sqw_args['Qmin'] = [2,1,1]
sqw_args['Qmax'] = [6,1,1]
params = Parameters.params(**md_args,**sqw_args)
calc = Calculator.calc(params)
io = FileIO.io(params,calc,f_name='sqe_7x7x7_5.dat')

sqw_args['Qmin'] = [2,0.25,0]
sqw_args['Qmax'] = [6,0.25,0]
params = Parameters.params(**md_args,**sqw_args)
calc = Calculator.calc(params)
io = FileIO.io(params,calc,f_name='sqe_7x7x7_6.dat')

sqw_args['Qmin'] = [2,0.25,0.25]
sqw_args['Qmax'] = [6,0.25,0.25]
params = Parameters.params(**md_args,**sqw_args)
calc = Calculator.calc(params)
io = FileIO.io(params,calc,f_name='sqe_7x7x7_7.dat')

