

import numpy as np


class c_rutile:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,basis=None,types=None,lattice_vectors=None):

        """
        define the primitive unitcell;

        nothing is really error checked here
        """

        # types of atoms
        if types is None:
            self.types = np.array([0,0,1,1,1,1])
        else:
            self.types = np.array(types)

        self.num_basis = self.types.size
    
        # positions of atoms in reduced coordinates
        if basis is None:
            self.basis = np.array([[0.5000000000000000,  0.5000000000000000,  0.5000000000000000],
                                   [0.0000000000000000,  0.0000000000000000,  0.0000000000000000],
                                   [0.1953400114833092,  0.8046599885166907,  0.5000000000000000],
                                   [0.8046599885166907,  0.1953400114833092,  0.5000000000000000],
                                   [0.3046599885166907,  0.3046599885166907,  0.0000000000000000],
                                   [0.6953400114833093,  0.6953400114833093,  0.0000000000000000]])
        else:
            self.basis = np.array(basis)

        # lattice vectors for crystal translations
        if lattice_vectors is None:
            self.lattice_vectors = np.array([[4.593,0.000,0.000],
                                             [0.000,4.593,0.000],
                                             [0.000,0.000,2.959]])
        else:
            self.lattice_vectors = np.array(lattice_vectors)

    # ----------------------------------------------------------------------------------------------

    def make_supercell(self,reps=[1,1,1]):

        """
        create rectangular supercell by replicating primitive cell reps[ii] number of times 
        along the ii^th direction 
        """

        self.reps = np.array(reps)
        self.num_reps = np.prod(self.reps)
        
    # ----------------------------------------------------------------------------------------------
