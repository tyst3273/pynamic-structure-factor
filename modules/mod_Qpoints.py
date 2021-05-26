import numpy as np
from mod_utils import print_stdout

class Qpoints:

    """
    store, intialize, and distribute Qpoints
    """

    # -----------------------------------------------------------------------------------------

    def generate_Qpoints(self,invars):
        """ 
        generate Q points, either a 2d scan or read list from file.
        """
        if invars.Qpoints_file == False:
            self._Qpoints_from_slice(invars) 
        else:
            self._Qpoints_from_list(invars)    

        message = f'Qpoints: {self.total_Qsteps}'
        print_stdout(message,msg_type='Brillouin zone path')

        for Q in range(self.total_Qsteps):
            message = (f'{Q+1}\t{self.total_reduced_Q[Q,0]:2.3f} {self.total_reduced_Q[Q,1]:2.3f} '
                    f'{self.total_reduced_Q[Q,2]:2.3f} r.l.u.')
            print_stdout(message)

    # --------------------------------------------------------------------------------------------

    def distribute_Q_over_procs(self,invars,num_ranks):
        """
        num_Q_per_proc is determined as the largest integer diving the total_Qsteps number. 
        the remainder is placed on rank 0 (if there is a remainder...)
        should probably use scatter/gather methods from MPI, but don't know how yet
        """
        # if only 1 proc, return the total Q set
        if num_ranks == 1:
            message = f'process: 0    Q points: {self.total_Qsteps:g}\n'
            print_stdout(message,msg_type='Q points on each process')
            self.Q_on_procs = [self.total_reduced_Q]
            return 

        # else, take num_Q_per_proc to be the largest integer dividing total_Qsteps
        num_Q_per_proc = self.total_Qsteps//num_ranks

        nQpp_arr = np.zeros(num_ranks).astype(int)
        nQpp_remainder = self.total_Qsteps-num_ranks*num_Q_per_proc

        nQpp_arr[:] = num_Q_per_proc
        nQpp_arr[0] = nQpp_arr[0]+nQpp_remainder

        if nQpp_remainder != 0:
            message = (f'one proc. will compute {nQpp_remainder+num_Q_per_proc} Qpoints;'
                        f' the rest will compute {num_Q_per_proc}.\n this might be inefficient.')
            print_stdout(message,msg_type='WARNING (PARALLELISM)')
        else:
            message = f'each proc. will compute {num_Q_per_proc} Qpoints'
            print_stdout(message,msg_type='NOTE (PARALLELISM)')

        message = f'process: 0    Q points: {nQpp_arr[0]:g}\n'
        for ii in range(1,num_ranks):
            message = message+f' process: {ii:g}    Q points: {nQpp_arr[ii]:g}\n'
        print_stdout(message,msg_type='Q points on each process')

        # put the Qpoints for each proc in the list
        self.Q_on_procs = []
        shift = 0
        for ii in range(num_ranks):
            self.Q_on_procs.append(self.total_reduced_Q[shift:shift+nQpp_arr[ii],:])
            shift = shift+nQpp_arr[ii]

    # ---------------------------------------------------------------------------------------------

    def rank_init(self,lattice,rank):
        """
        get the Qpoints for this rank and convert to 1/Angstrom
        """
        self.reduced_Q = np.array(self.Q_on_procs[rank])
        self.Qsteps = self.reduced_Q.shape[0]
        self._convert_Q_to_1_over_Angstrom(lattice)

    # ---------------------------------------------------------------------------------------------

    def reconvert_Q_points(self,lattice):
        self._convert_Q_to_1_over_Angstrom(lattice)

    # =======================================================================================
    # ------------------------------ private methods ----------------------------------------
    # =======================================================================================

    def _Qpoints_from_slice(self,invars):
        """
        generate Q-points between Qmin and Qmax for 2d BZ scan. 
        there will be total_Qsteps number of them. generate list in rlu.
        designed for orthogonal lattice vecs, but set up to be 
        easily exteneded to non-orthogonal case. 
        """
        self.total_Qsteps = invars.total_Qsteps
        self.total_reduced_Q = np.zeros((self.total_Qsteps,3))
        self.total_reduced_Q[:,0] = np.linspace(invars.Qmin[0],invars.Qmax[0],invars.total_Qsteps)
        self.total_reduced_Q[:,1] = np.linspace(invars.Qmin[1],invars.Qmax[1],invars.total_Qsteps)
        self.total_reduced_Q[:,2] = np.linspace(invars.Qmin[2],invars.Qmax[2],invars.total_Qsteps)

    # ---------------------------------------------------------------------------------------------

    def _Qpoints_from_list(self,invars):
        """
        Give a csv file of Qpoints. 1 per line, each coord seperated by spaces. Overwrites other 
        definitions for Q slices.
        """
        self.total_reduced_Q = np.loadtxt(invars.Qpoints_file)

        if len(self.total_reduced_Q.shape) == 1:
            self.total_reduced_Q = self.total_reduced_Q.reshape((1,3))

        self.total_Qsteps = self.total_reduced_Q.shape[0]

    # ---------------------------------------------------------------------------------------------

    def _convert_Q_to_1_over_Angstrom(self,lattice):
        """
        convert Q points from rlu to 1/A
        """
        self.Qpoints = np.zeros((self.Qsteps,3))
        for Q in range(self.Qsteps):
            self.Qpoints[Q,:] = (lattice.r_lattice_vectors[0,:]*self.reduced_Q[Q,0]+
                    lattice.r_lattice_vectors[1,:]*self.reduced_Q[Q,1]+
                    lattice.r_lattice_vectors[2,:]*self.reduced_Q[Q,2])







