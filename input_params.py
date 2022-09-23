
# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = \
    '/home/ty/research/projects/md_simulations/rutile/matsui_akaogi/sqw_xray/tmp/pos.h5'
# '/home/ty/research/projects/psf_data/pos.h5'
unwrap_trajectory = True
recalculate_box_lengths = False

# options for splitting up trajectory
num_trajectory_blocks = 5
trajectory_blocks = [0] #None 

# options for writing results
output_directory = None
output_prefix = 'psf'
calc_bragg = True
calc_diffuse = True
calc_sqw = True

# simulation options
md_time_step = 16 # femtoseconds; time step IN FILE, not used in simulation
md_num_steps = 625
md_num_atoms = 6480
md_supercell_reps = [6,6,30] 

# unit cell/ crystal info
lattice_vectors = [[4.593,0.000,0.000], # angstroms
                   [0.000,4.593,0.000],
                   [0.000,0.000,2.959]]
atom_types = ['Ti','Ti','O','O','O','O']
basis_positions = [[0.5000000000000000,  0.5000000000000000,  0.5000000000000000],
                   [0.0000000000000000,  0.0000000000000000,  0.0000000000000000],
                   [0.1953400114833092,  0.8046599885166907,  0.5000000000000000],
                   [0.8046599885166907,  0.1953400114833092,  0.5000000000000000],
                   [0.3046599885166907,  0.3046599885166907,  0.0000000000000000],
                   [0.6953400114833093,  0.6953400114833093,  0.0000000000000000]]

# experiment info
experiment_type = 'xrays' # 'neutrons' or 'xrays'

# options for how to generate Q-points for calculation
Qpoints_option = 'text_file' # mesh, mesh_file, write_mesh, text_file, or path

# for 'Qpoints_option' == 'mesh' ; 
# note, symmetry only works with plane centered on Q=(0,0,0) right now
# and requires spglib !!
Q_mesh_symmetry = False 
Q_mesh_H = [-2,2,40]
Q_mesh_K = [-2,2,40]
Q_mesh_L = 2

# 'Qpoints_option' == 'file'
Q_file = 'Qpts.dat'

# 'Qpoints_option' == 'path'
Q_path_start = [[0,0,0],
                [1,0,2]]
Q_path_end = [[1,0,2],
              [3,0,2]]
Q_path_steps = [20,20]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 1






