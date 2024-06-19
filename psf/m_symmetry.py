
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

        # call spglib depending on spin-pol
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
        msg += f'\nspacegroup found by spglib (w/o magnetism):\n  {self.spacegroup}'
        msg += f'\nnumber of symmetry operations (w/ magnetism):\n  {self.num_sym_ops}'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def get_irreducible_set_of_Qpoints(self):

        """
        loop over the set of Qpoints and find the irreducible set that is mapped onto all others
        
        let S be a rotation.

        since r' = S@r, k.r' = k.S@r = k@S.r, i.e. we look at kpts such that k' = k@S
        """

        _Qpts = self.comm.Qpoints
        _Q = _Qpts.Q_rlu
        _num_Q = _Qpts.num_Q

        for qq in range(_num_Q):
            
            print(qq)

        crash()

    # ----------------------------------------------------------------------------------------------

    
