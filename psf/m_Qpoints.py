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
from copy import deepcopy

from psf.m_error import crash

# --------------------------------------------------------------------------------------------------

class c_Qpoints:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm,timers):

        """
        holds lattice and reciprocal lattice vectors
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config
        self.timers = timers

        self.Qpoints_option = config.Qpoints_option
        
        # a convenience flag
        self.use_mesh = False
        if self.Qpoints_option == 'mesh':
            self.use_mesh = True

        self.use_symmetry = self.comm.symmetry.use_symmetry

    # ----------------------------------------------------------------------------------------------

    def generate_Qpoints(self):

        """
        sets Qpoints depending on what was set in the file

        ALL METHODS MUST PRODUCE 
            Q_rlu_full: np.array with shape [num_Q_full]x[3] and
            num_Q_full: int == number of Q-points
        """

        print('\n*** Q-points ***')

        if self.Qpoints_option == 'file':
            """
            read Q-points (in rlu) from txt file 
            """
            self._read_Q_text_file()

        elif self.Qpoints_option == 'path':
            """
            generate Q-points along paths thru reciprocal space. 
            """
            self._get_Q_along_paths()

        elif self.Qpoints_option == 'mesh':
            """
            generate a mesh in reciprocal space to calculate on
            """
            self._get_Q_on_mesh()
        
        # also need Qpts in cartesian coords
        self.get_cartesian_Qpoints()

        # now reduce the set using symmetry ...
        if self.use_symmetry:
            self.comm.symmetry.get_irreducible_set_of_Qpoints()
            self.Q_rlu = self.Q_rlu_full[self.irr_inds]
            self.Q_cart = self.Q_cart_full[self.irr_inds]
            self.num_Q = self.irr_inds.size
        else:
            self.Q_rlu = self.Q_rlu_full
            self.Q_cart = self.Q_cart_full
            self.num_Q = self.num_Q_full

        # print the Q-points to user
        msg = f'\nnumber of Q-points: {self.num_Q:g}\n'
        if self.num_Q >= 50:
            _loop = 50
            msg += 'only printing the first 50 Q-points!\n'
        else:
            _loop = self.num_Q

        msg += '          ------------- (rlu) -------------   ------------ (cart) -----------' 
        for ii in range(_loop):
            _ = f'\n Q[{ii:g}]:'
            msg += f'{_:10}'
            msg += f'{self.Q_rlu[ii,0]: 10.5f} {self.Q_rlu[ii,1]: 10.5f}' \
                    f' {self.Q_rlu[ii,2]: 10.5f}   ' 
            msg += f'{self.Q_cart[ii,0]: 10.5f} {self.Q_cart[ii,1]: 10.5f}' \
                    f'{self.Q_cart[ii,2]: 10.5f}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _get_Q_along_paths(self):

        """
        generate Q-points along path(s) thru BZ
        """

        msg = f'generating Q-points along paths'
        print(msg)

        _start = self.config.Q_path_start
        _end = self.config.Q_path_end
        _steps = self.config.Q_path_steps
        _num_paths = self.config.num_Q_paths 
        
        self.num_Q_full = np.sum(_steps)
        self.Q_rlu_full = np.zeros((self.num_Q_full,3),dtype=float)

        _shift = 0
        for ii in range(_num_paths):
            self.Q_rlu_full[_shift:_shift+_steps[ii],0] = np.linspace(_start[ii,0],
                                                    _end[ii,0],_steps[ii]+1)[:-1]
            self.Q_rlu_full[_shift:_shift+_steps[ii],1] = np.linspace(_start[ii,1],
                                                    _end[ii,1],_steps[ii]+1)[:-1]
            self.Q_rlu_full[_shift:_shift+_steps[ii],2] = np.linspace(_start[ii,2],
                                                    _end[ii,2],_steps[ii]+1)[:-1]
            _shift += _steps[ii]

    # ----------------------------------------------------------------------------------------------

    def get_cartesian_Qpoints(self):

        """
        get Qpoints in cartesian coords
        """

        self.Q_cart_full = self.comm.lattice.convert_coords(self.Q_rlu_full,
                self.comm.lattice.reciprocal_lattice_vectors)
        self.Q_len_full = np.sqrt(np.sum(self.Q_cart_full**2,axis=1))

    # ----------------------------------------------------------------------------------------------

    def _read_Q_text_file(self):

        """
        read Qpoints from text file    
        """

        self.Q_file = self.config.Q_file
        msg = f'getting Q-points from text file:\n  \'{self.Q_file}\''
        print(msg)

        try:
            _Q = np.loadtxt(self.Q_file,dtype=float)
        except Exception as _ex:
            msg = 'reading Q-points text file failed!\n'
            crash(msg,_ex)

        msg = 'Q-points in text file should have shape [num_Q]x[3]\n'
        _Q = np.atleast_2d(_Q)
        if _Q.shape[1] != 3:
            crash(msg)

        self.Q_rlu_full = _Q
        self.num_Q_full = self.Q_rlu_full.shape[0]

        msg = f'there were {self.num_Q_full:g} Q-points in the file'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _Q_steps(self,num_Q,Q_range):
        
        """
        return steps along Q axis based on mesh args
        """
        
        if num_Q == 1:
            Q = np.array([Q_range[0]])
        else:
            Q = np.linspace(Q_range[0],Q_range[1],num_Q) 

        return Q
    
    # ----------------------------------------------------------------------------------------------

    def _Q_range(self,Q_mesh):
        
        """
        parse Q_mesh_* arg and find number of Q-points and range. also return
        wheter it is symmetric about Q=0
        """
        
        _eps = 1e-6
        
        if Q_mesh.size == 1:
            num_Q = 1
            Q_range = np.array([Q_mesh,Q_mesh])
        else:
            num_Q = int(Q_mesh[2])
            Q_range = np.array([Q_mesh[0],Q_mesh[1]])

        return num_Q, Q_range

    # ----------------------------------------------------------------------------------------------

    def _get_Q_on_mesh(self):

        """
        generate Q-points on mesh that is specified by user
        """

        msg = 'generating Q-points on mesh'
        print(msg)

        self.Q_mesh_H = self.config.Q_mesh_H
        self.Q_mesh_K = self.config.Q_mesh_K
        self.Q_mesh_L = self.config.Q_mesh_L

        # get number and range of Q on each axis
        self.num_H, self.H_range = self._Q_range(self.Q_mesh_H)
        self.num_K, self.K_range = self._Q_range(self.Q_mesh_K)
        self.num_L, self.L_range = self._Q_range(self.Q_mesh_L)

        # for 'reshaping' the intensity grids
        self.mesh_shape = [self.num_H,self.num_K,self.num_L]

        # get steps along each Q axis
        self.H = self._Q_steps(self.num_H,self.H_range)
        self.K = self._Q_steps(self.num_K,self.K_range)
        self.L = self._Q_steps(self.num_L,self.L_range)

        self._get_full_Q_mesh()

    # ----------------------------------------------------------------------------------------------

    def _get_full_Q_mesh(self):

        """
        make Q-point mesh on full grid without using spglib
        """

        _H, _K, _L = np.meshgrid(self.H,self.K,self.L,indexing='ij')
        _H = _H.flatten(); _K = _K.flatten(); _L = _L.flatten()

        self.num_Q_full = _H.size
        self.Q_rlu_full = np.zeros((self.num_Q_full,3),dtype=float)
        self.Q_rlu_full[:,0] = _H; self.Q_rlu_full[:,1] = _K; self.Q_rlu_full[:,2] = _L

    # ----------------------------------------------------------------------------------------------

    def unfold_onto_cartesian_mesh(self,arr=None):

        """
        put data onto Q-points 'mesh'
        """
        
        _shape = list(arr.shape)
        if len(_shape) != 1:
            _mesh_shape = deepcopy(self.mesh_shape)
            _mesh_shape.extend(_shape[1:])
        else:
            _mesh_shape = deepcopy(self.mesh_shape)

        arr.shape = _mesh_shape
        
        return arr

    # ----------------------------------------------------------------------------------------------




