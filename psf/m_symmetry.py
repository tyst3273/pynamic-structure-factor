
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

        # fix me!
        _latvecs = self.config.symm_lattice_vectors
        _pos = self.config.symm_basis_pos
        _nums = self.config.symm_basis_types
        _magmoms = self.config.symm_magmoms

        cell = (_latvecs,_basis_pos,_basis_types)

        #rotations, translations, equivalent_atoms = spglib.get_symmetry(cell) 
        symm = spglib.get_symmetry(cell)
        print(symm)

    # ----------------------------------------------------------------------------------------------

    
