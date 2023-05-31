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


trajectory_format = 'lammps_hdf5' 

trajectory_file = 'pos.h5'

unwrap_trajectory = True

trajectory_stride = 1

trajectory_skip = 0

trajectory_trim = 0

num_trajectory_blocks = 1

trajectory_blocks = None

output_directory = None

output_prefix = 'psf'

calc_sqw = True

md_time_step = 16

md_num_steps = 6250

md_num_atoms = 6480

lattice_vectors = [[4.593,0.000,0.000], 
                   [0.000,4.593,0.000],
                   [0.000,0.000,2.959]]

box_vectors = None

atom_types = ['Ti','O']

experiment_type = 'neutrons' 

Qpoints_option = 'mesh' 

Q_mesh_H = [-3,3,41]

Q_mesh_K = [-3,3,41]

Q_mesh_L = 1

Q_file = 'Qpts.dat'

Q_path_start = [[0,0,0], 
                [2,0,0]]

Q_path_end = [[2,0,0],
              [2,0,2]]

Q_path_steps = [21,21]

num_Qpoint_procs = 1



