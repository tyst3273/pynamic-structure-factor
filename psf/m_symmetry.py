#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                                           !
#   ! Copyright 2025 by Tyler C. Sterling                                       !
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

from psf.m_error import crash
import numpy as np

# --------------------------------------------------------------------------------------------------

class c_symmetry:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm,timers):

        """
        use spglib to get space group, symmetry operations, etc. for given lattice vectors, 
        basis positions, and optionally spins
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config
        self.timers = timers

        self.symprec = 1e-5
        
        self.use_symmetry = self.config.use_symmetry
        if self.use_symmetry:
            self._get_symmetry()

    # ----------------------------------------------------------------------------------------------

    def _get_symmetry(self):

        """
        make calls to spglib. note, to enforce consistency with Q-points, the lattice vectors 
        are taken to be the ones used for Q-points. define the rest of the unitcell consistent
        with these
        """

        try:
            import spglib
        except Exception as _ex:
            crash(f'spglib couldnt be imported!',_ex)

        _symprec = self.symprec

        # need this stuff
        _latvecs = self.config.lattice_vectors
        _pos = self.config.symm_basis_pos
        _nums = self.config.symm_basis_types
        _magmoms = self.config.symm_magmoms

        # internation symbol / number
        self.spacegroup = spglib.get_spacegroup((_latvecs,_pos,_nums),symprec=_symprec)
        _version = spglib.spg_get_version()

        # call spglib depending on magnetism
        if _magmoms is None:
            cell = (_latvecs,_pos,_nums)
            symm = spglib.get_symmetry(cell=cell,symprec=_symprec)
        else:
            cell = (_latvecs,_pos,_nums,_magmoms)
            symm = spglib.get_magnetic_symmetry(cell=cell,symprec=_symprec)

        # these are in reduced coords
        self.rotations = symm['rotations']
        self.translations = symm['translations']

        self.equivalent_atoms = symm['equivalent_atoms']
        self.num_sym_ops = self.translations.shape[0]

        msg = '\n*** symmetry ***'
        msg += f'\nspglib: {_version}'
        msg += f'\nspacegroup found by spglib (w/o magnetism):\n  {self.spacegroup}'
        msg += f'\nnumber of symmetry operations (w/ magnetism):\n  {self.num_sym_ops}'
        print(msg)

        warn = '\n*** WARNING ***\n'
        warn += 'using symmetry isnt fully tested yet and I cant gaurantee there\n'
        warn += 'are no errors! please use caution and check results carefully!'
        print(warn)

    # ----------------------------------------------------------------------------------------------

    def get_irreducible_set_of_Qpoints(self):

        """
        loop over the set of Qpoints and find the irreducible set that is mapped onto all others
        
        let W be a rotation.
            since r' = W@r, k.r' = k.W@r = k@W.r, i.e. we look at kpts such that k' = k@W

        WRONG:
            since r' = W@r, k^T.r' = k^T.W@r = k^T@W.r, then k^T' = k^T@W, i.e. k' = W^T@k
            need to check if this matters!

        NOTE: this one is slower than the _get_irreducible... method below, but I am more confident
        that is correct

        NOTE: if this proves too slow, I could parallelize it over sym-ops

        NOTE: maybe I should only be using the symmomorphic ops? (i.e. no fractional tralslations)
        """

        self.timers.start_timer('get_irreducible_Q',units='s')

        msg = '\ngetting irreducible set of Q-points.\nthis might take a while ...'
        print(msg)

        # Q-points
        _Qpts = self.comm.Qpoints
        _Qrlu = _Qpts.Q_rlu_full
        _num_Q = _Qpts.num_Q_full

        # convert Qrlu to strings to make finding vectors in _Qrlu easier
        # maybe not the fastest way, but lets see if it works!
        _Qrlu_str = np.empty(_num_Q,dtype='<U23')
        for ii in range(_num_Q):
            _Qrlu_str[ii] = self._get_Q_str(_Qrlu[ii,:])

        # symmetry operations
        _num_sym = self.num_sym_ops
        _rot = self.rotations
        
        # this maps to irreducible set
        _map = -1*np.ones(_num_Q,dtype=int)

        for ii in range(_num_Q):

            if ii % 1000 == 0:
                print(f'  now on Qpt {ii}/{_num_Q}')
            
            # this Q is already mapped onto a previous one
            if _map[ii] >= 0:
                continue 

            for jj in range(_num_sym):  

                _W = _rot[jj,...] # rotation matrix
                _Qp = _Qrlu[ii,:]@_W # Q' = Q@W

                # find inds Q is mapped onto by this rotation
                _inds = self._find_Q_prime_in_Q_array(_Qp,_Qrlu_str)

                # this rotation doesnt map Q to any others
                if _inds is None:
                    continue

                # this rotation maps Q to Q[_ind] 
                for _ind in _inds:
                    _map[_ind] = ii
            
            # if map is complete, break
            if np.all(_map > 0):
                break

        # inds in the irreducible set, map from full to irreducible set
        _Qpts.irr_inds, _Qpts.irr_map = np.unique(_map,return_inverse=True)

        msg = f'\nnumber of Q-points in the full set: {_num_Q}'
        msg += f'\nnumber of Q-points in the irreducible set: {_Qpts.irr_inds.size}'
        print(msg)

        self.timers.stop_timer('get_irreducible_Q')
        
    # ----------------------------------------------------------------------------------------------

    def _get_irreducible_set_of_Qpoints(self):

        """
        loop over the set of Qpoints and find the irreducible set that is mapped onto all others
        
        let W be a rotation.
            since r' = W@r, k.r' = k.W@r = k@W.r, i.e. we look at kpts such that k' = k@W

        WARNING: this algorithm is faster, but I am not certain it is correct yet ...
        """

        self.timers.start_timer('get_irreducible_Q',units='s')

        msg = 'getting irreducible set of Q-points...'
        print(msg)

        # Q-points
        _Qpts = self.comm.Qpoints
        _Qrlu = _Qpts.Q_rlu
        _num_Q = _Qpts.num_Q

        # convert Qrlu to strings to make finding vectors in _Qrlu easier
        # maybe not the fastest way, but lets see if it works!
        _Qrlu_str = np.empty(_num_Q,dtype='<U23')
        for ii in range(_num_Q):
            _Qrlu_str[ii] = self._get_Q_str(_Qrlu[ii,:])

        # symmetry operations
        _num_sym = self.num_sym_ops
        _rot = self.rotations

        # this maps to irreducible set
        _map = -1*np.ones(_num_Q,dtype=int)

        # skip ii == 0, the idenity operation
        for ii in range(1,_num_sym):
            
            if ii % 10 == 0:
                print(f'  now on sym. op. {ii}/{_num_sym}')

            for jj in range(_num_Q):
                
                # if this Qpt is already mapped to another, skip
                if _map[jj] >= 0:
                    continue

                _W = _rot[ii,...]
                _Q = _Qrlu[jj,...]
                _Qp = _Q@_W 

                _inds = self._find_Q_prime_in_Q_array(_Qp,_Qrlu_str)

                if _inds is None:
                    continue

                for _ind in _inds:
                    if _map[_ind] >= 0:
                        continue
                    _map[_ind] = jj

            if np.all(_map > 0):
                break

        # ... now need to get unique array that maps to other, see get_irreducible... above

        self.timers.stop_timer('get_irreducible_Q')

    # ----------------------------------------------------------------------------------------------

    def _get_Q_str(self,Q):

        """
        convert Q 
        """
        
        return f'{Q[0]:.3f}_{Q[1]:.3f}_{Q[2]:.3f}'

    # ----------------------------------------------------------------------------------------------

    def _find_Q_prime_in_Q_array(self,Qp,Qrlu_str):
        
        """
        should allow finding multiple matching indices since it's allowed to calculate S(Q,w) on
        duplicate Q-points
        """

        _Qp_str = self._get_Q_str(Qp)
        _inds = np.flatnonzero(_Qp_str == Qrlu_str)

        if _inds.size == 0:
            return None

        return _inds

    # ----------------------------------------------------------------------------------------------

    def unfold_onto_full_Q_set(self,sq):

        """
        put data onto full set of Q-points
        """

        _Qpts = self.comm.Qpoints
        _num_Q_full = _Qpts.num_Q_full
        _irr_map = _Qpts.irr_map

        _ndim = sq.ndim

        # elastic has no energy dimension
        if _ndim == 1:
            sq_full = np.zeros((_num_Q_full),dtype=float)
            for ii in range(_num_Q_full):
                _ind = _irr_map[ii]
                sq_full[ii] = sq[_ind]
        
        # inelastic does
        else:
            sq_full = np.zeros((_num_Q_full,sq.shape[1]),dtype=float)
            for ii in range(_num_Q_full):
               _ind = _irr_map[ii]
               sq_full[ii,:] = sq[_ind,:]

        return sq_full

    # ----------------------------------------------------------------------------------------------

