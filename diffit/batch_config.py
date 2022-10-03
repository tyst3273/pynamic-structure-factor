
# options for where to get data/preprocessing
trajectory_format = 'external'
unwrap_trajectory = False

# options for splitting up trajectory
num_trajectory_blocks = 1
trajectory_blocks = None 

# options for writing results
output_directory = 'out' #None
output_prefix = 'psf'
calc_bragg = True
calc_rho_squared = False
calc_diffuse = True
calc_sqw = False

# simulation options
md_num_steps = 1

# unit cell/ crystal info
lattice_vectors = [[4.593,0.000,0.000], # angstroms
                   [0.000,4.593,0.000],
                   [0.000,0.000,2.959]]
atom_types = ['Ti','O']

# experiment info
experiment_type = 'neutrons' # 'neutrons' or 'xrays'

# options for how to generate Q-points for calculation
Qpoints_option = 'mesh' # mesh, mesh_file, write_mesh, text_file, or path

# for 'Qpoints_option' == 'mesh' ; 
Q_mesh_H = [0,4,100]
Q_mesh_K = [0,4,100]
Q_mesh_L = 0

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 16






