

# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5'
trajectory_file = '../lammps/pos.h5'

unwrap_trajectory = True

# options for splitting up trajectory
num_trajectory_blocks = 2
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
Qpoints_option = 'path' # mesh, file, or path

# 'Qpoints_option' == 'path'
Q_path_start = [[0,0,0],
                [5,0,0],
                [2.5,0,0]]
Q_path_end = [[5,0,0],
              [5,5,0],
              [2.5,5,0]]
Q_path_steps = [60,60,60]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 16






