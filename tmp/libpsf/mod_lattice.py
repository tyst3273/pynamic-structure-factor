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

import numpy as np
from libpsf.mod_utils import print_stdout

class lattice:

    """
    lattice and reciprocal lattice vectors 
    """

    # -----------------------------------------------------------------------------------------

    def __init__(self,invars):

        """
        store lattice and reciprocal lattice vectors
        """

        self.lattice_vectors = invars.lattice_vectors
        self.r_lattice_vectors = np.zeros((3,3))

        # print cell lengths from INPUT file
        message = (f'cell lengths from input: {self.lattice_vectors[0,0]} '
                   f'{self.lattice_vectors[1,1]} '
                   f'{self.lattice_vectors[2,2]} Angstrom')
        print_stdout(message,msg_type='NOTE')

        # print whether or not lattice vectors will be recalculated from traj. file
        if not invars.recalculate_cell_lengths:
            message = 'using cell lengths from input\n'
            print_stdout(message,msg_type='NOTE')
        else:
            message = 'using cell lengths from hdf5 trajectory file'
            print_stdout(message,msg_type='NOTE')

        # set up reciprocal lattice
        self._compute_reciprocal_lattice()

        # print the lattice/reciprocal lattice
        message = (f'real space lattice from input file (Angstrom):\n'
                f'  {self.lattice_vectors[0,0]: 2.3f} {self.lattice_vectors[0,1]: 2.3f}'
                f' {self.lattice_vectors[0,2]: 2.3f}\n  {self.lattice_vectors[1,0]: 2.3f}'
                f' {self.lattice_vectors[1,1]: 2.3f} {self.lattice_vectors[1,2]: 2.3f}\n'
                f'  {self.lattice_vectors[2,0]: 2.3f} {self.lattice_vectors[2,1]: 2.3f}'
                f' {self.lattice_vectors[2,2]: 2.3f}\n')
        print_stdout(message,msg_type='NOTE')

        message = (f'reciprocal space lattice from input file (1/Angstrom):\n'
                f'  {self.r_lattice_vectors[0,0]: 2.3f} {self.r_lattice_vectors[0,1]: 2.3f}'
                f' {self.r_lattice_vectors[0,2]: 2.3f}\n  {self.r_lattice_vectors[1,0]: 2.3f}'
                f' {self.r_lattice_vectors[1,1]: 2.3f} {self.r_lattice_vectors[1,2]: 2.3f}\n'
                f'  {self.r_lattice_vectors[2,0]: 2.3f} {self.r_lattice_vectors[2,1]: 2.3f}'
                f' {self.r_lattice_vectors[2,2]: 2.3f}')
        print_stdout(message)
            
    # ------------------------------------------------------------------------------------------

    def recompute_lattice(self):

        """
        recompute lattice vectors etc. from data read from traj file. 
        """

        self.r_lattice_vectors = np.zeros((3,3))
        self._compute_reciprocal_lattice()

        # print the lattice/reciprocal lattice
        message = (f'real space lattice from trajectory file (Angstrom):\n'
                f'  {self.lattice_vectors[0,0]: 2.3f} {self.lattice_vectors[0,1]: 2.3f}'
                f' {self.lattice_vectors[0,2]: 2.3f}\n  {self.lattice_vectors[1,0]: 2.3f}'
                f' {self.lattice_vectors[1,1]: 2.3f} {self.lattice_vectors[1,2]: 2.3f}\n'
                f'  {self.lattice_vectors[2,0]: 2.3f} {self.lattice_vectors[2,1]: 2.3f}'
                f' {self.lattice_vectors[2,2]: 2.3f}\n')
        print_stdout(message,msg_type='NOTE')

        message = (f'reciprocal space lattice from trajectory file (1/Angstrom):\n'
                f'  {self.r_lattice_vectors[0,0]: 2.3f} {self.r_lattice_vectors[0,1]: 2.3f}'
                f' {self.r_lattice_vectors[0,2]: 2.3f}\n  {self.r_lattice_vectors[1,0]: 2.3f}'
                f' {self.r_lattice_vectors[1,1]: 2.3f} {self.r_lattice_vectors[1,2]: 2.3f}\n'
                f'  {self.r_lattice_vectors[2,0]: 2.3f} {self.r_lattice_vectors[2,1]: 2.3f}'
                f' {self.r_lattice_vectors[2,2]: 2.3f}')
        print_stdout(message)

    # =======================================================================================
    # ------------------------------ private methods ----------------------------------------
    # =======================================================================================

    def _compute_reciprocal_lattice(self):

        """
        compute reciprocal lattice vectors from real lattice
        """

        self.cell_vol = self.lattice_vectors[0,:].dot(np.cross(self.lattice_vectors[1,:],
            self.lattice_vectors[2,:]))
        self.r_lattice_vectors[0,:] = 2*np.pi*np.cross(self.lattice_vectors[1,:],
                self.lattice_vectors[2,:])/self.cell_vol
        self.r_lattice_vectors[1,:] = 2*np.pi*np.cross(self.lattice_vectors[2,:],
                self.lattice_vectors[0,:])/self.cell_vol
        self.r_lattice_vectors[2,:] = 2*np.pi*np.cross(self.lattice_vectors[0,:],
                self.lattice_vectors[1,:])/self.cell_vol

    # ---------------------------------------------------------------------------------------












