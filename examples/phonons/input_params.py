
# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = '../lammps/pristine/pos.h5'

unwrap_trajectory = True

# options for splitting up trajectory
num_trajectory_blocks = 20
trajectory_blocks = [0] #None 

# options for writing results
output_directory = 'out' #None
output_prefix = 'psf'
calc_bragg = True
calc_rho_squared = True
calc_diffuse = True
calc_sqw = True

# simulation inputs
md_time_step = 24 # femtoseconds; time step IN FILE, not used in simulation
md_num_steps = 4000
md_num_atoms = 13824
md_supercell_reps = [12,12,12] 

# unit cell used to define Q-points in cartesian coords
lattice_vectors = [[ 5.431, 0.000, 0.000], # angstroms
                   [ 0.000, 5.431, 0.000],
                   [ 0.000, 0.000, 5.431]]

# vectors spanning the simulation cell
box_vectors = None

# experiment info
experiment_type = 'neutrons' # 'neutrons' or 'xrays'
atom_types = ['Si']

# options for how to generate Q-points for calculation
Qpoints_option = 'mesh' # mesh, file, or path

# for 'Qpoints_option' == 'mesh' ; 
Q_mesh_H = [-4,4,49]
Q_mesh_K = [-4,4,49]
Q_mesh_L = 0 #[-2,2,25]


# 'Qpoints_option' == 'file'
Q_file = 'Qpts.dat'

# 'Qpoints_option' == 'path'
Q_path_start = [[0,0,0],
                [1,0,2]]
Q_path_end = [[1,0,2],
              [3,0,2]]
Q_path_steps = [20,20]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 16





