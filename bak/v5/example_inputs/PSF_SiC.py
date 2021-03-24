import sys
sys.path.append('modules')

import Parameters
import Calculator
import Plot
import FileIO

# these are the parameters during the simulation
# dt is in seconds, b and lattice vectors in angstrom
md_args = {'traj_file':'lammps/SiC/NPT/pos_npt.dat',
           'file_format':'default',
           'dt':0.5e-15,
           'stride':16,
           'total_steps':2**18,
           'num_atoms':1280,
           'num_types':2,
           'supercell':[10,4,4],
           'b':{'1':4.1491*1e-5,'2':6.6460*1e-5},
           'lattice_vectors':[[4.326,0,0],
                              [0,4.326,0],
                              [0,0,4.326]]}

# these are parameters used to compute S(Q,w)
# Qmin, Qmax are in reciprocal lattice units
# Qstep is the number of steps along the BZ slice
sqw_args = {'Qmin':[2,0,0],
            'Qmax':[6,0,0],
            'Qsteps':80,
            'blocks':5,
            'debug':False}

params = Parameters.params(**md_args,**sqw_args)
calc = Calculator.calc(params)
#plot = Plot.color_map(params,calc,mask_tol=1e-8)
io = FileIO.io(params,calc,f_name='sqe_SiC_NPT.dat')

