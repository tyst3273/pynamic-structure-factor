
from psf.m_error import crash

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

        self._get_symmetry()

    # ----------------------------------------------------------------------------------------------

    def _get_symmetry(self):

        """
        make calls to spglib
        """

        try:
            import spglib
        except Exception as _ex:
            crash(f'spglib couldnt be imported!',_ex)

        _latvecs = self.comm.lattice.lattice_vectors
        _basis_pos = self.comm.atoms.atom_positions_reduced
        _basis_types = self.comm.atoms.atom_type_nums

        cell = (_latvecs,_basis_pos,_basis_types)

        #rotations, translations, equivalent_atoms = spglib.get_symmetry(cell) 
        symm = spglib.get_symmetry(cell)
        print(symm)

    # ----------------------------------------------------------------------------------------------

    
