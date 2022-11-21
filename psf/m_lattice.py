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

from psf.m_error import crash

class c_lattice:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm):

        """
        holds lattice and reciprocal lattice vectors
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config

        self.set_lattice_vectors(config.lattice_vectors)

    # ----------------------------------------------------------------------------------------------

    def set_lattice_vectors(self,lattice_vectors):

        """
        resets lattice vectors and recalculates reciprocal lattice vectors
        """

        self.lattice_vectors = lattice_vectors

        # check if numericall orthorhombic
        self._check_ortho()

        # get reciprocal lattice vectors
        self._get_reciprocal_lattice()

        msg = '\n*** lattice ***\n'
        msg += 'direct_lattice [Angstrom]'
        for ii in range(3):
            msg += '\n  '
            for jj in range(3):
                msg += f'{self.lattice_vectors[ii,jj]: 8.5f} '
        msg += '\nreciprocal_lattice [2*pi/Angstrom]'
        for ii in range(3):
            msg += '\n  '
            for jj in range(3):
                msg += f'{self.reciprocal_lattice_vectors[ii,jj]: 8.5f} '
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _check_ortho(self):

        """
        check if lattice vectors are numerically orthorhombic 
        """
       
        self.ortho_lattice_vectors = True
        for ii in range(3):
            for jj in range(3):
                if ii == jj:
                    continue
                else:
                    if np.round(self.lattice_vectors[ii,jj],6) != 0.0:
                        self.ortho_lattice_vectors = False
                        
        if not self.ortho_lattice_vectors:

            if self.config.unwrap_trajectory:
                msg = 'only orthogonal lattice vectors allowed in unwrap_trajectory=True\n'
                msg += f'if this functionality is really needed, contact the author'
                msg += ' at\n --- ty.sterling@colorado.edu\n'
                crash(msg)
            else:
                msg = 'non-orthogonal lattice vectors detected. this is still under development!\n'
                print('\n** warning **\n'+msg)

    # ----------------------------------------------------------------------------------------------

    def convert_coords(self,vecs,matrix):

        """
        take vectors and convert from one coordinate system to the other using the matrix given 
        as arg
        """

        _res = np.zeros(vecs.shape,dtype=float)
        _num = vecs.shape[0]
        
        for ii in range(_num):
            _res[ii,:] = vecs[ii,0]*matrix[0,:]+vecs[ii,1]*matrix[1,:]+vecs[ii,2]*matrix[2,:]

        return _res

    # ----------------------------------------------------------------------------------------------

    def _get_reciprocal_lattice(self):

        """
        get reciprocal lattice vectors
        """

        self.reciprocal_lattice_vectors = np.zeros((3,3),dtype=float)

        self.cell_vol = self.lattice_vectors[0,:].dot(np.cross(self.lattice_vectors[1,:],
            self.lattice_vectors[2,:]))
        self.reciprocal_lattice_vectors[0,:] = 2*np.pi*np.cross(self.lattice_vectors[1,:],
                self.lattice_vectors[2,:])/self.cell_vol
        self.reciprocal_lattice_vectors[1,:] = 2*np.pi*np.cross(self.lattice_vectors[2,:],
                self.lattice_vectors[0,:])/self.cell_vol
        self.reciprocal_lattice_vectors[2,:] = 2*np.pi*np.cross(self.lattice_vectors[0,:],
                self.lattice_vectors[1,:])/self.cell_vol

    # ----------------------------------------------------------------------------------------------


        
