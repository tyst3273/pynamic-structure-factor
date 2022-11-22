
# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = 'pos.h5'

# '/home/ty/research/projects/psf_data/pos.h5'
unwrap_trajectory = True

# options for splitting up trajectory
num_trajectory_blocks = 10
trajectory_blocks = None 

# options for writing results
output_directory = None
output_prefix = 'constQ'
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
experiment_type = 'neutrons' 

# options for how to generate Q-points for calculation
Qpoints_option = 'text_file' 

# 'Qpoints_option' == 'file'
Q_file = 'Qpts'

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 4






