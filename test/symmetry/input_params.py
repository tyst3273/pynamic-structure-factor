
"""
example input file
"""

# options for where to get data/preprocessing
trajectory_format = 'user_hdf5' 
trajectory_file = 'test.hdf5'
output_prefix = 'mesh'

unwrap_trajectory = False


# options for splitting up trajectory
num_trajectory_blocks = 5
trajectory_blocks = None 

# simulation options
md_time_step = 24 # femtoseconds; time step IN FILE, not used in simulation
md_num_steps = 2500
md_num_atoms = 13824

# unit cell/ crystal info
lattice_vectors = [[5.431,0.000,0.000], # angstroms
                   [0.000,5.431,0.000],
                   [0.000,0.000,5.431]]
atom_types = ['Si']

# experiment info
experiment_type = 'neutrons' # 'neutrons' or 'xrays'

# options for how to generate Q-points for calculation
Qpoints_option = 'mesh' # mesh, mesh_file, write_mesh, text_file, or path

# for 'Qpoints_option' == 'mesh' ; 
Q_mesh_H = [-6,6,145]
Q_mesh_K = [-6,6,145]
Q_mesh_L = 2

# 'Qpoints_option' == 'path'
Q_path_start = [[-1,1.25,0],
                [1,1.25,0],
                [1,-1.25,0]]
Q_path_end = [[1,1.25,0],
              [1,-1.25,0],
              [-1,-1.25,0]]
Q_path_steps = [24,18,24]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 24

# symmetry
use_symmetry = True
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



