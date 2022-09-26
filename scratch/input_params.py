
# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = '/home/ty/research/projects/md_simulations/rutile/matsui_akaogi/sqw_xray/pos.h5' 
#trajectory_format = 'external'
unwrap_trajectory = True

# options for splitting up trajectory
num_trajectory_blocks = 5
trajectory_blocks = [0] #None 

# options for writing results
output_directory = 'out' #None
output_prefix = 'psf'
calc_bragg = True
calc_rho_squared = True
calc_diffuse = True
calc_sqw = True

# simulation options
md_time_step = 1 # femtoseconds; time step IN FILE, not used in simulation
md_num_steps = 100000//16
md_num_atoms = 6*6*6*30
md_supercell_reps = [6,6,30] 

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
Q_mesh_H = [0,4,40]
Q_mesh_K = [0,4,40]
Q_mesh_L = 1

# 'Qpoints_option' == 'file'
Q_file = 'Qpts.dat'

# 'Qpoints_option' == 'path'
Q_path_start = [[2.5,2.5,0.5],
                [2.5,2.5,3.5],
                [3.5,3.5,3.5]]
Q_path_end = [[2.5,2.5,3.5],
              [3.5,3.5,3.5],
              [0,0,0]]
Q_path_steps = [90,12,36]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 16






