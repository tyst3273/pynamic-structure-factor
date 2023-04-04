

# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = '../lammps/pristine/pos.h5'

unwrap_trajectory = True

# options for splitting up trajectory
num_trajectory_blocks = 4
trajectory_blocks = 0

# options for writing results
output_directory = 'out' #None
output_prefix = 'pristine'
calc_bragg = True
calc_rho_squared = True
calc_diffuse = True
calc_sqw = True

# simulation inputs
md_time_step = 16     # femtoseconds; time step IN FILE, not Verlet time step
md_num_steps = 2000   # number of steps IN FILE
md_num_atoms = 13824 

#md_num_atoms = 13074 # number of atoms in vacancy cell

# unit cell used to define Q-points in cartesian coords
lattice_vectors = [[ 4.593, 0.000, 0.000], # angstroms
                   [ 0.000, 4.593, 0.000],
                   [ 0.000, 0.000, 2.959]]

# vectors spanning the simulation cell
box_vectors = None

# experiment info
experiment_type = 'neutrons' # 'neutrons' or 'xrays'
atom_types = ['Ti','O']

# options for how to generate Q-points for calculation
Qpoints_option = 'mesh' # mesh, file, or path

# for 'Qpoints_option' == 'mesh' ; 
Q_mesh_H = [0,6,73]
Q_mesh_K = [0,6,73]
Q_mesh_L = 2


# number of processes to split Q-point parallelization over
num_Qpoint_procs = 16






