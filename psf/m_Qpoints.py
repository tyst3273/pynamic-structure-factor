
import numpy as np

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

    # ----------------------------------------------------------------------------------------------

    def generate_Qpoints(self):

        """
        sets Qpoints depending on what was set in the file

        ALL METHODS MUST PRODUCE 
            Q_rlu: np.array with shape [num_Q]x[3] and
            num_Q: int == number of Q-points

        """

        self.timers.start_timer('generate_Qpoints',units='ms')

        if self.Qpoints_option == 'text_file':

            """
            read Q-points (in rlu) from txt file 
            """
            self.Q_file = self.config.Q_file
            self._read_Q_text_file()

        elif self.Qpoints_option == 'mesh_file':

            """
            read Q-points and other mesh data from hdf5 file.
            can be used to calculate on same mesh as another calculation
            """
            self.Q_file = self.config.Q_file
            #self._read_Q_mesh_file()
            crash('\'mesh_file\' is not available yet!\n')

        elif self.Qpoints_option == 'path':
            
            """
            generate Q-points along paths thru reciprocal space. 
            """
            self._get_Q_along_paths()

        elif self.Qpoints_option == 'mesh' or self.Qpoints_option == 'write_mesh':

            """
            generate a mesh in reciprocal space to calculate on. can use symmetry to reduce
            number of points, but only if the mesh is centered on Q=(0,0,0). otherwise, 
            arbitrary mesh can be done without using symmetry
            """
            self._get_Q_on_mesh()

            if self.Qpoints_option == 'write_mesh':
                self.Q_file = self.config.Q_file
                #self._write_mesh_file()
                crash('\'write_mesh\' is not available yet!\n')
        
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

        self.timers.stop_timer('generate_Qpoints')

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

    def _Q_steps(self,Q_mesh,eps=0.0001):
        
        """
        return steps along Q axis based on mesh args
        """
        
        if Q_mesh.size == 1:
            num_Q = 1
            Q = np.copy(Q_mesh).reshape(1,)
            _0 = True
        else:
            num_Q = int(Q_mesh[2])
            Q = np.linspace(Q_mesh[0],Q_mesh[1],num_Q)
            if np.abs(Q.mean()) < eps:
                _0 = True
            else:
                _0 = False

        return num_Q, Q, _0

    # ----------------------------------------------------------------------------------------------

    def _get_Q_on_mesh(self):

        """
        generate Q-points on mesh that is specified by user

        note that grid has to be centered on Q=(0,0,0) for symmetry reduction to work. actually,
        if cutting 2D planes, the axes spanning the plane have to be centered on 0... i think
        it's fine to ignore the center of the otherone, since all of that Q_i in the set are the 
        same. 
        """

        msg = '\ngenerating Q-points on mesh'
        print(msg)

        self.Q_mesh_symmetry  = self.config.Q_mesh_symmetry
        self.Q_mesh_H = self.config.Q_mesh_H
        self.Q_mesh_K = self.config.Q_mesh_K
        self.Q_mesh_L = self.config.Q_mesh_L

        # get steps along each Q axis
        self.num_H, self.H, _H0 = self._Q_steps(self.Q_mesh_H)
        self.num_K, self.K, _K0 = self._Q_steps(self.Q_mesh_K)
        self.num_L, self.L, _L0 = self._Q_steps(self.Q_mesh_L)

        # grid has to be centered on Q=(0,0,0) for symmetry reduction to work 
        if _H0 and _K0 and _L0:
            _0_centered = True
        else:
            _0_centered = False
        
        # try to grid using symmetry to reduce number of points
        if self.Q_mesh_symmetry:
            if not _0_centered:
                msg = 'to use Q-point symmetry, grid must be centered on Q=(0,0,0)\n' \
                      'either disable symmetry or pick a different grid!\n'
                crash(msg)
            if self.config.basis_positions is None:
                msg = 'basis_positions must be defined to find space group and\n' \
                       'reduce number of Q-points using symmetry\n'
                crash(msg)
            self._get_symmetry_reduced_Q_mesh()

        # otherwise use full grid
        else:
            self._make_Q_mesh()
    
    # ----------------------------------------------------------------------------------------------

    def _make_Q_mesh(self):

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

    def _get_symmetry_reduced_Q_mesh(self):

        """
        generate Q-points on mesh using spglib to reduce number of points
        """

        try:
            import spglib
        except Exception as _ex:
            msg = 'spglib couldnt be imported. check the installation or disable symmetry\n' \
                  'for reducing the number of Q-points. see the exception below for details.\n'
            crash(msg,_ex)

        lat_vecs = self.comm.lattice.lattice_vectors
        basis_pos = self.config.basis_positions
        basis_num = self.config.type_map
        cell = (lat_vecs,basis_pos,basis_num)

        mesh = [self.num_H,self.num_K,self.num_L]
        mapping, grid = spglib.get_ir_reciprocal_mesh(mesh,cell,is_shift=[0,0,0])

        # full grid 
        self.num_Q_full = mapping.size
        self.Q_rlu_full = grid/np.array(mesh,dtype=float)

        # put on extended Q-grid, not just in 1st BZ
        self.Q_rlu_full[:,0] *= 2*self.H.max()
        self.Q_rlu_full[:,1] *= 2*self.K.max()
        self.Q_rlu_full[:,2] *= 2*self.L.max()
        
        # get irreducible part
        irr_inds, self.Q_wts = np.unique(mapping,return_counts=True)
        self.Q_wts = self.Q_wts/self.num_Q_full
        self.Q_rlu = self.Q_rlu_full[irr_inds]
        self.num_Q = self.Q_rlu.shape[0]

        # map from irreducible part to full grid
        self.Q_irr_to_full = np.zeros(mapping.shape,dtype=int)
        for ii in range(self.num_Q_full):
            self.Q_irr_to_full[ii] == np.flatnonzero(mapping[ii] == irr_inds)[0]

        msg = 'number of Q-points on full grid:\n'
        msg += f'  {self.num_Q_full:g}\n'
        msg += 'numer of Q-point in irreducible part:\n'
        msg += f'  {self.num_Q:g}\n'
        print(msg)

    # ----------------------------------------------------------------------------------------------



