
# how much info to print
verbosity = 1

# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = 'pos.h5'
unwrap_trajectory = True
recalculate_box_lengths = False

# options for writing results
output_directory = None
output_prefix = 'psf'
calc_bragg = True
calc_timeavg = True
calc_sqw = True

# simulation options
md_time_step = 16 # femtoseconds
md_num_steps = 6250
md_num_atoms = 6480
md_supercell_reps = [6,6,30] 

# unit cell/ crystal info
lattice_vectors = [[4.593,0.000,0.000], # angstroms
                   [0.000,4.593,0.000],
                   [0.000,0.000,2.959]]
atom_types = ['Ti','Ti','O','O','O','O']
basis_positions = None


# experiment info
experiment_type = 'neutrons' # 'neutrons' or 'xrays'


# options for how to generate Q-points for calculation
Qpoints_option = 'file' # mesh, file, or path

# 'Qpoints_option' == 'mesh'
Q_mesh_symmetry = True # used for 'mesh' option; requires spglib !!
Q_mesh = [48,48,1] 
Q_mesh_H = [-3,3]
Q_mesh_K = [-3,3]
Q_mesh_L = 1

# 'Qpoints_option' == 'file'
Q_file = 'Qpts.dat'

# 'Qpoints_option' == 'path'
Q_path_start = [[0,0,0], 
               [2,0,0]]
Q_path_end = [[2,0,0],
             [2,0,2]]
Q_path_steps = [21,21]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 1






