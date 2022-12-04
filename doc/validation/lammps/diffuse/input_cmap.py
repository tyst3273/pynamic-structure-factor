
# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = '../pos.h5'

# '/home/ty/research/projects/psf_data/pos.h5'
unwrap_trajectory = True

# options for splitting up trajectory
num_trajectory_blocks = 10
trajectory_blocks = None 

# options for writing results
output_directory = None
output_prefix = 'diffuse'
calc_bragg = True
calc_rho_squared = True
calc_diffuse = True
calc_sqw = True

# simulation options
md_time_step = 24 
md_num_steps = 100000//24
md_num_atoms = 12**3*8
md_supercell_reps = [12,12,12] 

# unit cell/ crystal info
lattice_vectors = [[5.431,0.000,0.000], 
                   [0.000,5.431,0.000],
                   [0.000,0.000,5.431]]
atom_types = ['Si']

# experiment info
experiment_type = 'xrays' 

# options for how to generate Q-points for calculation
Qpoints_option = 'mesh'
Q_mesh_H = [-4,4,12*8]
Q_mesh_K = [-4,4,12*8]
Q_mesh_L = 1 #[-4,4,12*8]

# symmetry
use_Qpoints_symmetry = True
symm_types = [1,1,1,1,1,1,1,1]
symm_lattice_vectors = [[1,0,0],
                        [0,1,0],
                        [0,0,1]]
symm_positions = [[0.00,0.00,0.00],
                  [0.00,0.50,0.50],
                  [0.50,0.00,0.50],
                  [0.50,0.50,0.00],
                  [0.25,0.25,0.25],
                  [0.75,0.75,0.25],
                  [0.25,0.75,0.75],
                  [0.75,0.25,0.75]]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 4






