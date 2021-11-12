
# system modules

# custom modules

class c_sqw:

    """
    class that does the calculating
    it calculates bragg, diffuse, and dynamic scattering if requested.
    """
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,logger,lattice,qpts,traj,xlens,timers):

        """
        set everything up. get frequencies, Q-spacing, etc. also 'allocate'
        sqw, bragg, and diffuse arrays
        """
        
        self.config = config
        self.logger = logger
        self.lattice = lattice
        self.qpts = qpts
        self.traj = traj
        self.xlens = xlens
        self.timers = timers

    # ----------------------------------------------------------------------------------------------

    def run_calc(self):

        """
        loop over blocks. for each, parse trajectories, unwrap if requested, then split the
        Q-points over num_process processes which run independently. once done, they write 
        the results to the Q-points in the bragg, sqw, diffuse arrays and merge back to calling
        process
        """

        msg = 'run_calc() not implemented yet'
        self.logger.warn(msg,trace=True)

    # ----------------------------------------------------------------------------------------------

