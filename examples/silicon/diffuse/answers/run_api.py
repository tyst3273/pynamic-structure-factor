
from psf.m_PSF import c_PSF


# common configuration options for all calculations are read from input file
input_file = 'input_pristine.py'
PSF = c_PSF(input_file)


# options can be over written here and passed as args below
num_trajectory_blocks = 2
trajectory_blocks = None
num_Qpoint_procs = 8


# --------------------------------------------------------------------------------------------------

# pristine neutrons
PSF.setup_calculation(trajectory_file = '../lammps/pristine/pos.h5',
                      output_prefix = 'pristine_neutrons',
                      atom_types = ['Si'],
                      experiment_type = 'neutrons',
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------

# pristine xrays
PSF.setup_calculation(trajectory_file = '../lammps/pristine/pos.h5',
                      output_prefix = 'pristine_xrays',
                      atom_types = ['Si'],
                      experiment_type = 'xrays',
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------

# substitutions neutrons
PSF.setup_calculation(trajectory_file = '../lammps/substitutions/pos.h5',
                      output_prefix = 'substitutions_neutrons',
                      atom_types = ['Si','Ge'],
                      experiment_type = 'neutrons',
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------

# substitutions xrays
PSF.setup_calculation(trajectory_file = '../lammps/substitutions/pos.h5',
                      output_prefix = 'substitutions_xrays',
                      atom_types = ['Si','Ge'],
                      experiment_type = 'xrays',
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------





