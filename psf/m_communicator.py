

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
        self.lattice = m_lattice.c_lattice(self.config,self.comm)

        # Q-points created using one of numerous methods specified in config
        self.Qpoints = m_Qpoints.c_Qpoints(self.config,self.comm,self.timers)
        self.Qpoints.generate_Qpoints()

        # scattering lengths for the S(Q,w) calculation
        self.xlengths = m_scattering_lengths.c_scattering_lengths(self.config,self.comm)


    # ----------------------------------------------------------------------------------------------


