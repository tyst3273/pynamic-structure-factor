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

# --------------------------------------------------------------------------------------------------

"""
description: to do ...
required: 
type:
"""

trajectory_format = 'lammps_hdf5' 

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

trajectory_file = 'pos.h5'

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

unwrap_trajectory = False

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

trajectory_stride = 1

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

trajectory_skip = 0

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

trajectory_trim = 0

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

num_trajectory_blocks = 1

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

trajectory_blocks = None

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

output_directory = None

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

output_prefix = 'psf'

# --------------------------------------------------------------------------------------------------


"""
description:
required: True
type:
"""

md_time_step = None # must be defined

# --------------------------------------------------------------------------------------------------

"""
description:
required: True
type:
"""

md_num_steps = None 

# --------------------------------------------------------------------------------------------------

"""
description:
required: True
type:
"""

md_num_atoms = None 

# --------------------------------------------------------------------------------------------------

"""
description:
required: True
type:
"""

lattice_vectors = None 

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

box_vectors = None

# --------------------------------------------------------------------------------------------------

"""
description:
required:  True
type:
"""

atom_types = None

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

experiment_type = 'neutrons' 

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

Qpoints_option = 'mesh' 

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

Q_mesh_H = [-3,3,41]

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

Q_mesh_K = [-3,3,41]

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

Q_mesh_L = 0

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

Q_file = 'Qpts.dat'

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

Q_path_start = [[0,0,0], 
                [2,0,0]]

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

Q_path_end = [[2,0,0],
              [2,0,2]]

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

Q_path_steps = [21,21]

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

num_Qpoint_procs = 1

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

use_symmetry = False

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

symm_basis_pos = None

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

symm_basis_types = None 

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

symm_magmoms = None

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

calc_incoherent = True

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

calc_coherent = True

# --------------------------------------------------------------------------------------------------

"""
description:
required: 
type:
"""

calc_sqw = True

# --------------------------------------------------------------------------------------------------









