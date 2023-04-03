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
description: format of the trajectory data file. 
    current supported args:  
    'lammps_hdf5'
    must be written using commands like:
       dump            POSITIONS all h5md 5 pos.h5 position species
       dump_modify     POSITIONS sort id
    note that lammps wraps the positions by default, so 'unwrap_trajectory' should
    be set to True
    'user_hdf5'
    is a standard format made by converting text files written by lammps,
    vasp, etc. the 'merge_trj.py' script can be used to convert the text files. 
     it automatically unwraps the trajectories so dont do it again.
type: str
"""
trajectory_format = 'lammps_hdf5' 

# --------------------------------------------------------------------------------------------------
"""
description: path to the trajectory file
type: str
"""
trajectory_file = 'pos.h5'

# --------------------------------------------------------------------------------------------------
"""
description: whether or not to 'unwrap' the trajectory. by default, lammps moves atoms
    that cross the simulation box boundary back to the other side. this can cause the the 
    trajectory to jump discontinuously by a box vector, causing issues with the debye waller
    factor and other things
type: bool
"""
unwrap_trajectory = True

# --------------------------------------------------------------------------------------------------
"""
description: only use every n^th step. useful if the sampling in the file is finer than needed
type: int
"""
trajectory_stride = 1

# --------------------------------------------------------------------------------------------------
"""
description: skip this many steps at the beginning of the trajectory. useful if 
    e.g. volume is equilibrated at the beginning of the file, so skip all the steps
    where fluctuations are large
type: int
"""
trajectory_skip = 0

# --------------------------------------------------------------------------------------------------
"""
description: skip this many steps at the end of the trajectory. 
type: int
"""
trajectory_trim = 0

# --------------------------------------------------------------------------------------------------
"""
description: split up the remaining trajectory (after skipping, trimming, and down sampling)
    into this many 'blocks'. S(Q,w) is calculated for each block and then averaged over the 
    the blocks given in 'trajectory_blocks' below. 
type: int
"""
num_trajectory_blocks = 10

# --------------------------------------------------------------------------------------------------
"""
description: which blocks to actually use. split trajectory in 'num_trajectory_blocks' number
    of blocks and calculate on the ones given in this list. e.g. only use every other block
    to make sure that the S(Q,w) for each block are uncorrelated
    if trajectory_blocks == None, all are used
type: list
"""
trajectory_blocks = None

# --------------------------------------------------------------------------------------------------
"""
description: where to write the results of the calculation. 
type: str
"""
output_directory = None

# --------------------------------------------------------------------------------------------------
"""
description: 'prefix' to prepend to the output file. 
type: str
"""
output_prefix = 'psf'

# --------------------------------------------------------------------------------------------------
"""
description: calculate bragg scattering: |<exp(iQ.r(t))>|**2
type: bool
"""
calc_bragg = True

# --------------------------------------------------------------------------------------------------
"""
description: calculate 'rho' squared: rho = exp(iQ.r(t)), rho squared = |exp(iQ.r(t))|**2
type: bool
"""
calc_rho_squared = True

# --------------------------------------------------------------------------------------------------
"""
description: calculate diffuse scattering: <|exp(iQ.r(t))|>**2
type: bool
"""
calc_diffuse = True

# --------------------------------------------------------------------------------------------------
"""
description: calculate S(Q,w): S(Q,w) =  | \int dt exp(iQ.r(t)-iwt) |**2
type: bool
"""
calc_sqw = True

# --------------------------------------------------------------------------------------------------
"""
description: time step in the file in femtoseconds.
    not the timestep in the simulation!
    and NOT accounting for 'trajectory_stride' above. 
    the effective step is md_time_step*trajectory_stride, but this is calculated inside 
    the code and the user doesnt have to care about it
type: int
"""
md_time_step = 16 

# --------------------------------------------------------------------------------------------------
"""
description: number of steps in the trajectory file
type: int
"""
md_num_steps = 6250

# --------------------------------------------------------------------------------------------------
"""
description: number of atoms in the simulation
type: int
"""
md_num_atoms = 6480

# --------------------------------------------------------------------------------------------------
"""
description: number or replications of unitcell in the supercell
type: list with 3 elements, each an integer
"""
md_supercell_reps = [6,6,30] 

# --------------------------------------------------------------------------------------------------
"""
description: lattice vectors of unitcell. these used to calculate the Q points and to 'unwrap'
    the trajectories. they can be arbitrary, but keep in mind the Q-points have to be commensurate
    with the supercell to get sensible data. note that these should be in whatever length units
    are in the trajectory file.
        these are different from 'prim_lattice_vectors' in that these are used to calculate Q 
    in cartesian coords; 'prim_lattice_vectors' are the idealized lattice vectors used to 
    determine the space group and reduce the Q-point set using symmetry
type: [3]x[3] list of floats
"""
lattice_vectors = [[4.593,0.000,0.000], 
                   [0.000,4.593,0.000],
                   [0.000,0.000,2.959]]

# --------------------------------------------------------------------------------------------------
"""
description: the atoms types in the simulation. the code expects the types are consecutive
    integers starting at 1 and that there is one for each type. note, these are used to set the
    scattering lengths and form factors. e.g. if you want to calculate scattering from 
    deuterium rather than hydrogen, set the type to 2H rather than H etc
type: list of str.
"""
atom_types = ['Ti','O']

# --------------------------------------------------------------------------------------------------
"""
description: the experiment type, either neutrons or x-rays.
type: str
"""
experiment_type = 'neutrons' 

# --------------------------------------------------------------------------------------------------
"""
description: how to generate the set of Q-points. currently supported args are:
    'mesh': calculate S(Q,w) on automatically generated 3D Q-point mesh. see Q_mesh_H, Q_mesh_K,
        and Q_mesh_L
    'path': calculate S(Q,w) on automatically generated path thru Q-space. see Q_path_start,
        Q_path_end, and Q_path_steps below
    'file': calculate S(Q,w) on set of Q-points in csv file. see Q_file 
type: str
"""
Qpoints_option = 'mesh' 

# --------------------------------------------------------------------------------------------------
"""
description: set up Q-points with H-component spanning from Q_mesh_*[0] to Q_mesh_*[1] with 
    Q_mesh_*[2] steps along this edge of the grid. for only one Q-point along this edge, just 
    give a list with a single integer in it.
type: list with either 3 integer elements or 1 element
"""
Q_mesh_H = [-3,3,41]

# --------------------------------------------------------------------------------------------------
"""
description: set up Q-points with H-component spanning from Q_mesh_*[0] to Q_mesh_*[1] with 
    Q_mesh_*[2] steps along this edge of the grid. for only one Q-point along this edge, just 
    give a list with a single integer in it.
type: list with either 3 integer elements or 1 element
"""
Q_mesh_K = [-3,3,41]

# --------------------------------------------------------------------------------------------------
"""
description: set up Q-points with H-component spanning from Q_mesh_*[0] to Q_mesh_*[1] with 
    Q_mesh_*[2] steps along this edge of the grid. for only one Q-point along this edge, just 
    give a list with a single integer in it.
type: list with either 3 integer elements or 1 element
"""
Q_mesh_L = 1

# --------------------------------------------------------------------------------------------------
"""
description: csv file with Q-points in them. it should have shape Nx3. note it is read by 
    numpy.loadtxt() so make sure it will work with that.
type: str
"""
Q_file = 'Qpts.dat'

# --------------------------------------------------------------------------------------------------
"""
description: beginning of 'vertices' on Q-path.
type: Nx3 lists.
"""
Q_path_start = [[0,0,0], 
               [-2,0,0]]

# --------------------------------------------------------------------------------------------------
"""
description: end of vertices on Q-path. if you want consecutive paths, then Q_path_end[i] should
    be that same as Q_path_start[i-1]. should have same shape as Q_path_start
type: Nx3 list
"""
Q_path_end = [[2,0,0],
             [25,0,-2]]

# --------------------------------------------------------------------------------------------------
"""
description: number of steps along each segment of the Q-path. should have same number of 
    segments as Q_path_start and Q_path_end
type: N element list of integers
"""
Q_path_steps = [21,21]

# --------------------------------------------------------------------------------------------------
"""
description: number of processes to split Q-points over. Q-points are split-up round-robin
    style and then processes are spawned using multiprocessing library. note that python doesnt
    really support shared memory, so use this with caution if the trajectory is large (in the
    sense of how much disk space it requires)
type: int
"""
num_Qpoint_procs = 1

# --------------------------------------------------------------------------------------------------
"""
description: used to determine symmetry and reduce the Q-point. they are integer types 
    corresponding to the atoms in 'symm_positions'. default is unset
type: list of ints
example:
    symm_atom_types = [0,0,1,1,1,1] # 'Ti', 'Ti', 'O', 'O', 'O', 'O'
"""
symm_types = None 

# --------------------------------------------------------------------------------------------------
"""
description: reduced coordinates of atoms in unitcell corresponding to 'symm_lattice_vectors' 
    variable and to 'symm_types'. these are passed (with 'symm_lattice_vectors' and 
    'symm_types') to spglib to determine spacegroup. default is unset
type: nx3 list of floats.
example:
    symm_positions = [[0.5000000000000000,  0.5000000000000000,  0.5000000000000000],
                      [0.0000000000000000,  0.0000000000000000,  0.0000000000000000],
                      [0.1953400114833092,  0.8046599885166907,  0.5000000000000000],
                      [0.8046599885166907,  0.1953400114833092,  0.5000000000000000],
                      [0.3046599885166907,  0.3046599885166907,  0.0000000000000000],
                      [0.6953400114833093,  0.6953400114833093,  0.0000000000000000]]
"""
symm_positions = None 

# --------------------------------------------------------------------------------------------------
"""
description: idealized lattice vectors used to determine spacegroup. passed to spglib
    these are different from 'lattice_vectors' in that these are only used to determine
    that space group to reduce the Q-point set. 'lattice_vectors' is used to convert 
    Q to cartesian coords in 1/Angstrom. default is unset
type: [3]x[3] list of floats
example:
    symm_lattice_vectors = [[4.593,0.000,0.000], 
                            [0.000,4.593,0.000],
                            [0.000,0.000,2.959]]
"""
symm_lattice_vectors = None

# --------------------------------------------------------------------------------------------------
"""
description: use symmetry to reduce number of Q-points. must have mesh centered on Q=0 and must
    have spglib installed
type: bool
"""
use_Qpoints_symmetry = False







