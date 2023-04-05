

# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = '../lammps/pos.h5'

unwrap_trajectory = True

# options for splitting up trajectory
num_trajectory_blocks = 1
trajectory_blocks = None

# options for writing results
output_directory = 'out' #None
output_prefix = 'protonated'

# simulation inputs
md_time_step = 50     # femtoseconds; time step IN FILE, not Verlet time step
md_num_steps = 1000   # number of steps IN FILE
md_num_atoms = 20736 

#md_num_atoms = 13074 # number of atoms in vacancy cell

# unit cell used to define Q-points in cartesian coords
lattice_vectors = [[ 6.309, 0.000, 0.000], # angstroms
                   [ 0.000, 6.309, 0.000],
                   [ 0.000, 0.000, 6.309]]

# vectors spanning the simulation cell
box_vectors = None

# experiment info
experiment_type = 'neutrons' # 'neutrons' or 'xrays'
atom_types = ['C','N','H','H','Pb','I']

# options for how to generate Q-points for calculation
Qpoints_option = 'mesh' #'mesh' # mesh, file, or path

# for 'Qpoints_option' == 'mesh' ; 
Q_mesh_H = [-4,4,97]
Q_mesh_K = [-4,4,97]
Q_mesh_L = 2.5

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 16






