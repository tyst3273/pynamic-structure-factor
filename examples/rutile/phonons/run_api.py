
from PSF import c_PSF


# common configuration options for all calculations are read from input file
input_file = 'input_pristine.py'
PSF = c_PSF(input_file)


# options can be over written here and passed as args below
num_trajectory_blocks = 2
trajectory_blocks = None
num_Qpoint_procs = 16


# --------------------------------------------------------------------------------------------------

# pristine neutrons
PSF.setup_calculation(trajectory_file = '../lammps/pristine/pos.h5',
                      output_prefix = 'pristine_neutrons',
                      experiment_type = 'neutrons',
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------

# pristine xrays
PSF.setup_calculation(trajectory_file = '../lammps/pristine/pos.h5',
                      output_prefix = 'pristine_xrays',
                      experiment_type = 'xrays',
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------

vacancy_lattice_vectors = [[ 4.597, 0.000, 0.000], # angstroms
                           [ 0.000, 4.597, 0.000],
                           [ 0.000, 0.000, 2.961]]

# vacancies neutrons
PSF.setup_calculation(trajectory_file = '../lammps/vacancies/pos.h5',
                      output_prefix = 'vacancies_neutrons',
                      experiment_type = 'neutrons',
                      md_num_atoms = 13074,
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------

# vacancies xrays
PSF.setup_calculation(trajectory_file = '../lammps/vacancies/pos.h5',
                      output_prefix = 'vacancies_xrays',
                      experiment_type = 'xrays',
                      md_num_atoms = 13074,
                      num_trajectory_blocks = num_trajectory_blocks,
                      trajectory_blocks = trajectory_blocks,
                      num_Qpoint_procs = num_Qpoint_procs)
PSF.run()

# --------------------------------------------------------------------------------------------------





