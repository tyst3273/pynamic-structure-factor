
import numpy as np

from psf.m_error import crash

class c_Qpoints:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config):

        """
        holds lattice and reciprocal lattice vectors
        """

        self.Qpoints_option = config.Qpoints_option

    # ----------------------------------------------------------------------------------------------

    def set_Qpoints(self,config):

        """
        sets Qpoints depending on what was set in the file
        """

        if self.Qpoints_option == 'text_file':
            pass
        elif self.Qpoints_option == 'mesh_file':
            # read mesh from a file and calculate on mesh, unfold to full BZ, etc
            pass
        elif self.Qpoints_option == 'path':
            pass
        elif self.Qpoints_option == 'mesh':
            pass
        elif self.Qpoints_option == 'write_mesh':
            # write to Q-points mesh file and exit. used to chain together calculations
            pass

    # ----------------------------------------------------------------------------------------------


