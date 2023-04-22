
from psf.m_PSF import c_PSF


# common configuration options for all calculations are read from input file
input_file = 'input.py'
PSF = c_PSF(input_file)


# options can be over written here and passed as args below
num_trajectory_blocks = 2
trajectory_blocks = None
num_Qpoint_procs = 16



# --------------------------------------------------------------------------------------------------

# protonated neutrons
PSF.setup_calculation(trajectory_file = '../lammps/pos.h5',
                      output_prefix = 'protonated_neutrons',
                      experiment_type = 'neutrons',
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------

# deuterated neutrons
PSF.setup_calculation(trajectory_file = '../lammps/pos.h5',
                      output_prefix = 'deuterated_neutrons',
                      experiment_type = 'neutrons',
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs,
                      atom_types = ['C','N','2H','2H','Pb','I'])
PSF.run()

# --------------------------------------------------------------------------------------------------

# protonated xrays
PSF.setup_calculation(trajectory_file = '../lammps/pos.h5',
                      output_prefix = 'protonated_xrays',
                      experiment_type = 'xrays',
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------
