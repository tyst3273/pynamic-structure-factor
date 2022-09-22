
import psf.m_lattice as m_lattice
import psf.m_Qpoints as m_Qpoints
import psf.m_timing as m_timing
import psf.m_scattering_lengths as m_scattering_lengths
import psf.m_trajectory as m_trajectory
import psf.m_structure_factors as m_structure_factors
import psf.m_parallel as m_parallel


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
        self.Qpoints = m_Qpoints.c_Qpoints(self.config,self,self.timers)
        self.Qpoints.generate_Qpoints()

        # set up parallelism over Q-points
        self.paral = m_parallel.c_multi_processing(self.config,self)

        # find trajectory file and set up for block-averaging
        self.traj = m_trajectory.c_trajectory(self.config,self,self.timers)

        # scattering lengths for the S(Q,w) calculation
        self.xlengths = m_scattering_lengths.c_scattering_lengths(self.config,self)

        # setup object to hold structure factors, calculate stuff, etc.
        self.strufacs = m_structure_factors.c_structure_factors(self.config,self,self.timers)

    # ----------------------------------------------------------------------------------------------


