
import numpy as np

from psf.m_error import crash

class c_lattice:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config):

        """
        holds lattice and reciprocal lattice vectors
        """

        self.ortho_lattice_vectors = config.ortho_lattice_vectors
        
        if config.recalculate_box_lengths:
            msg = '\nlattice vectors will be dynamically recalculate using data\n' \
                  'from the trajectory file!'
            print(msg)

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
        
        for ii in range(3):
            for jj in range(3):
                if ii == jj:
                    continue
                else:
                    if np.round(self.lattice_vectors[ii,jj],6) != 0.0:
                        self.ortho_lattice_vectors = False
                        
        if not self.ortho_lattice_vectors:
            msg = 'only orthogonal (i.e. diagonal) lattice_vectors allowed for now\n'
            msg += f'if this functionality is really needed, contact the author'
            msg += ' at\n --- ty.sterling@colorado.edu\n'
            crash(msg)

    # ----------------------------------------------------------------------------------------------

    def convert_coords(self,vecs,matrix):

        """
        take vectors and convert from one coordinate system to the other using the matrix given 
        as arg
        """

        pass

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


        
