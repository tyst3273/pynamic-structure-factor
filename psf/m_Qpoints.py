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
        self.use_Qpoints_symmetry = config.use_Qpoints_symmetry
        
        # a convenience flag
        self.use_mesh = False
        if self.Qpoints_option == 'mesh':
            self.use_mesh = True

    # ----------------------------------------------------------------------------------------------

    def generate_Qpoints(self):

        """
        sets Qpoints depending on what was set in the file

        ALL METHODS MUST PRODUCE 
            Q_rlu: np.array with shape [num_Q]x[3] and
            num_Q: int == number of Q-points

        """

        if self.Qpoints_option == 'file':
            """
            read Q-points (in rlu) from txt file 
            """
            self.Q_file = self.config.Q_file
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

        # print the Q-points to user
        msg = f'\n*** Q-points ***\nnumber of Q-points: {self.num_Q:g}\n'
        if self.num_Q >= 50:
            loop = 50
            msg += 'only printing the first 50!\n'
        else:
            loop = self.num_Q

        msg += '          ---------- (rlu) ---------   --------- (cart) --------' 
        for ii in range(loop):
            _ = f'\n Q[{ii:g}]:'
            msg += f'{_:10}'
            msg += f'{self.Q_rlu[ii,0]: 8.5f} {self.Q_rlu[ii,1]: 8.5f}' \
                    f' {self.Q_rlu[ii,2]: 8.5f}   ' 
            msg += f'{self.Q_cart[ii,0]: 8.5f} {self.Q_cart[ii,1]: 8.5f}' \
                    f'{self.Q_cart[ii,2]: 8.5f}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _get_Q_along_paths(self):

        """
        generate Q-points along path(s) thru BZ
        """

        _start = self.config.Q_path_start
        _end = self.config.Q_path_end
        _steps = self.config.Q_path_steps
        _num_paths = self.config.num_Q_paths 
        
        self.num_Q = np.sum(_steps)
        self.Q_rlu = np.zeros((self.num_Q,3),dtype=float)

        _shift = 0
        for ii in range(_num_paths):
            self.Q_rlu[_shift:_shift+_steps[ii],0] = np.linspace(_start[ii,0],
                                                    _end[ii,0],_steps[ii]+1)[:-1]
            self.Q_rlu[_shift:_shift+_steps[ii],1] = np.linspace(_start[ii,1],
                                                    _end[ii,1],_steps[ii]+1)[:-1]
            self.Q_rlu[_shift:_shift+_steps[ii],2] = np.linspace(_start[ii,2],
                                                    _end[ii,2],_steps[ii]+1)[:-1]
            _shift += _steps[ii]

    # ----------------------------------------------------------------------------------------------

    def get_cartesian_Qpoints(self):

        """
        get Qpoints in cartesian coords
        """

        self.Q_cart = self.comm.lattice.convert_coords(self.Q_rlu,
                self.comm.lattice.reciprocal_lattice_vectors)
        self.Q_len = np.sqrt(np.sum(self.Q_cart**2,axis=1))

    # ----------------------------------------------------------------------------------------------

    def _read_Q_text_file(self):

        """
        read Qpoints from text file    
        """

        msg = f'\ngetting Q-points from text file:\n  \'{self.Q_file}\''
        print(msg)

        try:
            _Q = np.loadtxt(self.Q_file,dtype=float)
        except:
            msg = 'reading Q-points text file failed!\n'
            crash(msg)

        msg = 'Q-points in text file should have shape [num_Q]x[3]\n'
        if len(_Q.shape) == 1:
            if _Q.size != 3:
                crash(msg)
            else:
                _Q.shape = [1,3]

        self.Q_rlu = _Q
        self.num_Q = self.Q_rlu.shape[0]

        msg = f'there were {self.num_Q:g} Q-points in the file'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _Q_steps(self,num_Q,Q_range):
        
        """
        return steps along Q axis based on mesh args

        NOTE: consistent with what phonopy returns, this does NOT include the lowest
        Q point on the requested mesh
        """
        
        if num_Q == 1:
            Q = np.array([Q_range[0]])
        else:
            Q = np.linspace(Q_range[0],Q_range[1],num_Q+1)[1:] # trim off lowest point in mesh ...

        return Q
    
    # ----------------------------------------------------------------------------------------------

    def _Q_range(self,Q_mesh):
        
        """
        parse Q_mesh_* arg and find number of Q-points and range. also return
        wheter it is symmetric about Q=0
        """
        
        _eps = 1e-6
        
        symmetric = False
        
        if Q_mesh.size == 1:
            num_Q = 1
            Q_range = np.array([Q_mesh,Q_mesh])
            symmetric = True
        else:
            num_Q = int(Q_mesh[2])
            Q_range = np.array([Q_mesh[0],Q_mesh[1]])
            if np.abs(Q_range.mean()) < _eps:
                symmetric = True

        return num_Q, Q_range, symmetric

    # ----------------------------------------------------------------------------------------------

    def _get_Q_on_mesh(self):

        """
        generate Q-points on mesh that is specified by user
        """

        msg = '\ngenerating Q-points on mesh'
        print(msg)

        self.Q_mesh_H = self.config.Q_mesh_H
        self.Q_mesh_K = self.config.Q_mesh_K
        self.Q_mesh_L = self.config.Q_mesh_L

        # get number and range of Q on each axis
        self.num_H, self.H_range, self.H_symmetric = self._Q_range(self.Q_mesh_H)
        self.num_K, self.K_range, self.K_symmetric = self._Q_range(self.Q_mesh_K)
        self.num_L, self.L_range, self.L_symmetric = self._Q_range(self.Q_mesh_L)
        self.symmetric_mesh = bool(self.H_symmetric*self.K_symmetric*self.L_symmetric)
        
        # get steps along each Q axis
        self.H = self._Q_steps(self.num_H,self.H_range)
        self.K = self._Q_steps(self.num_K,self.K_range)
        self.L = self._Q_steps(self.num_L,self.L_range)

        if self.use_Qpoints_symmetry:
            self._get_Q_mesh_from_spglib()
        else:
            self._get_full_Q_mesh()

        # for 'reshaping' the intensity grids
        self.mesh_shape = [self.num_H,self.num_K,self.num_L]

    # ----------------------------------------------------------------------------------------------

    def _get_Q_mesh_from_spglib(self):

        """
        get Q-point mesh on reduced grid using spglib

        the Q-point mesh has to be symmetric about Q=0 to work since spglib creates 
        symmetric mesh. could probably work around that too but this is good enough for now
        
        NOTE: it is assumed that if an axis has only one Q point, it is symmetric about that 
        axis. this is unlikely to be correct... but is a hack for now. e.g. to calculate
        S(Q) in a plane perpendicular to L, the Q-point grid is reduced with L=0 and then 
        shifted such that L = whatever the user requested. 
        ... i should undo this in a later version
        """
        
        self.timers.start_timer('Q_from_spglib',units='s')

        try:
            import spglib
        except Exception as _ex:
            msg = 'spglib could not be imported! check the installation and try again\n' \
                  'or set use_Qpoints_symmetry = False. Printing the exception below.\n'
            crash(msg,exception=_ex)
            
        if not self.symmetric_mesh:
            msg = 'to use spglib to reduce number of Q-points, the mesh must be symmetric\n' \
                  'about Q=0. pick a new symmetric mesh of set use_Qpoints_symmetry=False\n'
            crash(msg)
            
        # spglib needs a tuple of lists
        _vecs = self.config.symm_lattice_vectors.tolist()
        _pos = self.config.symm_positions.tolist()
        _nums = self.config.symm_types
        _cell = (_vecs,_pos,_nums)
        _mesh = [self.num_H,self.num_K,self.num_L]

        # print spacegroup info for diagnostrics
        _spg = spglib.get_spacegroup(_cell,symprec=1e-5,angle_tolerance=-1.0,symbol_type=0)
        
        # call spglib to get irreducible mesh
        mapping, self.Q_rlu_full = spglib.get_ir_reciprocal_mesh(
                        _mesh,_cell,is_shift=[0,0,0])
        
        # full Qpoints mesh in rlu
        self.Q_rlu_full = self.Q_rlu_full.astype(float)
        self.Q_rlu_full[:,0] *= self.H_range[1]*2/_mesh[0]
        self.Q_rlu_full[:,1] *= self.K_range[1]*2/_mesh[1]
        self.Q_rlu_full[:,2] *= self.L_range[1]*2/_mesh[2]
        self.num_Q_full = self.Q_rlu_full.shape[0]

        # if the axes have only 1 Qpt, then we set that axis to be equal to the requested Q
        if self.num_H == 1:
            self.Q_rlu_full[:,0] = self.H[0]
        if self.num_K == 1:
            self.Q_rlu_full[:,1] = self.K[0]
        if self.num_L == 1:
            self.Q_rlu_full[:,2] = self.L[0]
        
        # inds of irr Q in full mesh and inds mapping full mesh to irr Q
        self.Q_irr_inds, self.Q_full_to_irr_map = np.unique(mapping,return_inverse=True)
        
        # irreducible Q-points in rlu
        self.num_Q = self.Q_irr_inds.size
        self.Q_rlu = self.Q_rlu_full[self.Q_irr_inds]

        msg = f'\nspace group found by spglib: {_spg}\n'
        msg += f'number of Q-points on full mesh: {self.num_Q_full}\n'
        msg += f'number of Q-points on irreducible mesh: {self.num_Q}'
        print(msg)
    
        self.timers.stop_timer('Q_from_spglib')
    
    # ----------------------------------------------------------------------------------------------

    def _get_full_Q_mesh(self):

        """
        make Q-point mesh on full grid without using spglib
        """

        _H, _K, _L = np.meshgrid(self.H,self.K,self.L,indexing='ij')
        _H = _H.flatten(); _K = _K.flatten(); _L = _L.flatten()

        self.num_Q = _H.size
        self.Q_rlu = np.zeros((self.num_Q,3),dtype=float)
        self.Q_rlu[:,0] = _H; self.Q_rlu[:,1] = _K; self.Q_rlu[:,2] = _L

        msg = 'number of Q-points on full grid:\n'
        msg += f'  {self.num_Q:g}\n'
        print(msg)
        
    # ----------------------------------------------------------------------------------------------
    
    def _put_on_full_mesh(self,arr):
        
        """
        unfold from reduced mesh on to full mesh
        """
        
        _shape = arr.shape
        if len(_shape) != 1:
            _full_arr = arr[self.Q_full_to_irr_map,:]
        else:
            _full_arr = arr[self.Q_full_to_irr_map]
        return _full_arr

    # ----------------------------------------------------------------------------------------------

    def unfold_onto_Q_mesh(self,arr=None):

        """
        put data onto Q-point on 'mesh'
        """
        
        if self.use_Qpoints_symmetry:
            arr = self._put_on_full_mesh(arr)

        _shape = list(arr.shape)
        if len(_shape) != 1:
            _mesh_shape = deepcopy(self.mesh_shape)
            _mesh_shape.extend(_shape[1:])
        else:
            _mesh_shape = deepcopy(self.mesh_shape)

        arr.shape = _mesh_shape
        
        if self.use_Qpoints_symmetry:
            arr = np.fft.fftshift(arr,axes=(0,1,2))

        return arr

    # ----------------------------------------------------------------------------------------------




