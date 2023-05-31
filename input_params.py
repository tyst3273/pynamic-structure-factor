
# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = '/home/ty/research/projects/md_simulations/rutile/matsui_akaogi/sqw_xray/pos.h5'

unwrap_trajectory = True

# options for splitting up trajectory
num_trajectory_blocks = 10
trajectory_blocks = None 

# options for writing results
output_directory = 'out' 
output_prefix = 'psf'

# simulation options
md_time_step = 16 # femtoseconds; time step IN FILE, not used in simulation
md_num_steps = 6250
md_num_atoms = 6480

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
Q_mesh_H = [-1,1,4]
Q_mesh_K = [-1,1,4]
Q_mesh_L = [-1,1,4]

# 'Qpoints_option' == 'file'
Q_file = 'Qpts.dat'

# 'Qpoints_option' == 'path'
Q_path_start = [[0,0,0],
                [1,0,2]]
Q_path_end = [[1,0,2],
              [3,0,2]]
Q_path_steps = [20,20]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 7






