#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                                           !
#   ! Copyright 2021 by Tyler C. Sterling and Dmitry Reznik,                    !
#   ! University of Colorado Boulder                                            !
#   !                                                                           !
#   ! This file is part of the pynamic-structure-factor (PSF) software.         !
#   ! PSF is free software: you can redistribute it and/or modify it under      !
#   ! the terms of the GNU General Public License as published by the           !
#   ! Free software Foundation, either version 3 of the License, or             !
#   ! (at your option) any later version. PSF is distributed in the hope        !
#   ! that it will be useful, but WITHOUT ANY WARRANTY; without even the        !
#   ! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  !
#   ! See the GNU General Public License for more details.                      !
#   !                                                                           !
#   ! A copy of the GNU General Public License should be available              !
#   ! alongside this source in a file named gpl-3.0.txt. If not see             !
#   ! <http://www.gnu.org/licenses/>.                                           !
#   !                                                                           !
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# options for where to get data/preprocessing
trajectory_format = 'lammps_hdf5' 
trajectory_file = \
        '/home/ty/research/projects/md_simulations/rutile/matsui_akaogi/sqw_xray/tmp/pos.h5'
unwrap_trajectory = True
recalculate_box_lengths = False

# options for splitting up trajectory
num_trajectory_blocks = 10
trajectory_blocks = [0,1,2,3,4,5,6,7,8,9]

# options for writing results
output_directory = None
output_prefix = 'psf'
calc_bragg = True
calc_rho_squared = True
calc_diffuse = True
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
Qpoints_option = 'mesh' # mesh, mesh_file, write_mesh, text_file, or path

# for 'Qpoints_option' == 'mesh' ;
# note, symmetry only works with plane centered on Q=(0,0,0) right now
# and requires spglib !!
Q_mesh_symmetry = True # used for 'mesh' option; requires spglib !!
Q_mesh_H = [-3,3,41]
Q_mesh_K = [-3,3,41]
Q_mesh_L = 1

# 'Qpoints_option' == 'file'
Q_file = 'Qpts.dat'

# 'Qpoints_option' == 'path'
Q_path_start = [[0,0,0], 
               [-2,0,0]]
Q_path_end = [[2,0,0],
             [25,0,-2]]
Q_path_steps = [21,21]

# number of processes to split Q-point parallelization over
num_Qpoint_procs = 1






