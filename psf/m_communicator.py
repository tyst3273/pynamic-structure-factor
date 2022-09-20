
import psf.m_lattice as m_lattice
import psf.m_Qpoints as m_Qpoints
import psf.m_timing as m_timing
import psf.m_scattering_lengths as m_scattering_lengths


class c_communicator:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,timers):
        
        """
        put other objects inside of this to conveniently pass them in/out of classes and methods
        """

        self.config = config
        self.timers = timers

    # ----------------------------------------------------------------------------------------------

    def setup_calculation(self):

        """
        initialize errything
        """

        # lattice and reciprocal lattice vectors
        self.lattice = m_lattice.c_lattice(self.config,self)

        # Q-points created using one of numerous methods specified in config
        self.Qpoints = m_Qpoints.c_Qpoints(self.config,self.timers,self)
        self.Qpoints.generate_Qpoints()

        # scattering lengths for the S(Q,w) calculation
        self.xlengths = m_scattering_lengths.c_scattering_lengths(self.config,self)

    # ----------------------------------------------------------------------------------------------


