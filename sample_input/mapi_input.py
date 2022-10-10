
# options for where to get data/preprocessing
trajectory_format = 'user_hdf5' 
trajectory_file = 'rappe_mapi_350K.hdf5'

# '/home/ty/research/projects/psf_data/pos.h5'
unwrap_trajectory = False

# options for splitting up trajectory
num_trajectory_blocks = 20
trajectory_blocks = [0,8,16] 

# options for writing results
output_directory = 'out' #None
output_prefix = 'psf'
calc_bragg = True
calc_rho_squared = True
calc_diffuse = True
calc_sqw = True

# simulation options
md_time_step = 50 # femtoseconds; time step IN FILE, not used in simulation
md_num_steps = 2001
md_num_atoms = 49152
md_supercell_reps = [16,16,16] 

# unit cell/ crystal info
lattice_vectors = [[6.295625,0.000000,0.000000], # angstroms
                   [0.000000,6.295625,0.000000],
                   [0.000000,0.000000,6.295625]]
atom_types = ['Pb','I','H','N','C']

# experiment info
experiment_type = 'xrays' # 'neutrons' or 'xrays'

# options for how to generate Q-points for calculation
Qpoints_option = 'mesh' # mesh, mesh_file, write_mesh, text_file, or path

# for 'Qpoints_option' == 'mesh' ; 
Q_mesh_H = [0,4,16*4]
Q_mesh_K = [0,4,16*4]
Q_mesh_L = 1.5

# 'Qpoints_option' == 'file'
Q_file = 'Qpts.dat'

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 4






