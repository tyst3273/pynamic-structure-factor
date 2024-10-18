
"""
example input file
"""

# options for where to get data/preprocessing
trajectory_format = 'user_hdf5' 
trajectory_file = 'pos.hdf5'
output_prefix = 'mesh'

unwrap_trajectory = False

# options for splitting up trajectory
num_trajectory_blocks = 8
trajectory_blocks = None 

# simulation options
md_time_step = 16 # femtoseconds; time step IN FILE, not used in simulation
md_num_steps = 4000
md_num_atoms = 7200
calc_sqw = True
calc_coherent = True
calc_incoherent = True 

# unit cell/ crystal info
lattice_vectors = [[4.5011,0.0000,0.0000], # angstroms
                   [0.0000,4.5011,0.0000],
                   [0.0000,0.0000,3.0182]]
atom_types = ['Ti','O']

# experiment info
experiment_type = 'neutrons' # 'neutrons' or 'xrays'

# options for how to generate Q-points for calculation
Qpoints_option = 'path' # mesh, mesh_file, write_mesh, text_file, or path

# for 'Qpoints_option' == 'mesh' ; 
Q_mesh_H = [3,4,13]
Q_mesh_K = [3,4,13]
Q_mesh_L = 2

# 'Qpoints_option' == 'path'
Q_path_start = [[0,0,0],
                [4,4,0],
                [8,0,0],
                [4,0,0],
                [4,0,4],
                [4,4,4]]
Q_path_end = [[4,4,0],
              [8,0,0],
              [4,0,0],
              [4,0,4],
              [4,4,4],
              [0,0,4]]
Q_path_steps = [40,40,40,48,40,40]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 12

# symmetry
use_symmetry = False
symm_basis_pos = [[0.00, 0.00, 0.00],
                  [0.00, 0.50, 0.50],
                  [0.50, 0.00, 0.50],
                  [0.50, 0.50, 0.00],
                  [0.25, 0.25, 0.25],
                  [0.75, 0.75, 0.25],
                  [0.75, 0.25, 0.75],
                  [0.25, 0.75, 0.75]]
symm_basis_types = [1,1,1,1,1,1,1,1]
#symm_magmoms = [0,0,0,0,0,0,0,0]



