#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                               !
#   ! this file is part of the 'pynamic-structure-factor' code      !
#   ! written by Ty Sterling at the University of Colorado          !
#   ! Boulder, advised by Dmitry Reznik.                            !
#   !                                                               !
#   ! the software calculates inelastic neutron dynamic structure   !
#   ! factors from molecular dynamics trajectories.                 !
#   !                                                               !
#   ! this is free software distrubuted under the GNU GPL v3 and    !
#   ! with no warrantee or garauntee of the results. you should     !
#   ! have recieved a copy of the new license with this software    !
#   ! if you do find bugs or have questions, dont hesitate to       !
#   ! write to the author at ty.sterling@colorado.edu               !
#   !                                                               !
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import numpy as np
from mod_utils import print_stdout, PSF_exception

class Qpoints:

    """
    store, intialize, and distribute Qpoints over procs
    """

    # -----------------------------------------------------------------------------------------

    def generate_Qpoints(self,invars):
        """ 
        generate Q points, either a 2d scan or read list from file.
        """
        # generate list of Q in rlu 
        if invars.Qpoints_file == False:
            self._Qpoints_from_slice(invars) # 2d scan from input file
        else:
            self._Qpoints_from_list(invars) # read Q points from file

        # print list of Q in rlu 
        message = f'Qpoints: {self.total_Qsteps}'
        print_stdout(message,msg_type='Brillouin zone path')
        Q_count = 0
        for Q in range(self.total_Qsteps):
            message = (f'{Q+1}\t{self.total_reduced_Q[Q,0]: 2.3f} {self.total_reduced_Q[Q,1]: 2.3f} '
                    f'{self.total_reduced_Q[Q,2]: 2.3f} r.l.u.')
            print_stdout(message)

            if Q_count == 49:
                message = ('...............................\n number of Qpoints is >= 50.'
                           '\n suppressing output.')
                print_stdout(message)
                break
            Q_count = Q_count+1

    # --------------------------------------------------------------------------------------------

    def distribute_Q_over_procs(self,invars,num_ranks):
        """
        num_Q_per_proc is determined as the largest integer dividing the total_Qsteps number. 
        the remainder is placed on rank 0 (if there is a remainder...)
        should probably use scatter/gather methods from MPI, but don't know how yet
        """

        # set up array of how many Q on each proc
        nQpp_arr = np.zeros(num_ranks).astype(int) # Num_Q_Per_Processor_ARRay
        proc = 0
        for qq in range(self.total_Qsteps):
            if proc == num_ranks:
                proc = 0
            nQpp_arr[proc] = nQpp_arr[proc]+1
            proc = proc+1

        # if any procs have 0 Q, print a warning
        if nQpp_arr.min() == 0:
            message = 'atleast one processor will have 0 Qpoints.\n this is probably not effecient.'
            print_stdout(message,msg_type='WARNING')

        # print parellelism info
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
        self.reduced_Q = np.array(self.Q_on_procs[rank]) # get the Q points for this rank
        self.Qsteps = self.reduced_Q.shape[0]            # get the number on this rank
        self._convert_Q_to_1_over_Angstrom(lattice)      # convert the Q from rlu to 1/A

    # ---------------------------------------------------------------------------------------------

    def reconvert_Q_points(self,lattice):
        """
        interface to use if recalcuting cell lenghts from traj file: the recip. lattice changes
        so we have to get Q in 1/A in the new basis
        """
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
        definitions for Q slices if a file name is given.
        """
        try:
            self.total_reduced_Q = np.loadtxt(invars.Qpoints_file) # read the Q points
        except:
            message = f'Qpoints file \'{invars.Qpoints_file}\' is broken'
            raise PSF_exception(message)

        if len(self.total_reduced_Q.shape) == 1: # if only 1 Q, reshape to avoid breaking stuff later
            self.total_reduced_Q = self.total_reduced_Q.reshape((1,3))

        self.total_Qsteps = self.total_reduced_Q.shape[0] # number of Q points

    # ---------------------------------------------------------------------------------------------

    def _convert_Q_to_1_over_Angstrom(self,lattice):
        """
        convert Q points from rlu to 1/A in cartesian coords
        """
        self.Qpoints = np.zeros((self.Qsteps,3))
        for Q in range(self.Qsteps):
            self.Qpoints[Q,:] = (lattice.r_lattice_vectors[0,:]*self.reduced_Q[Q,0]+
                    lattice.r_lattice_vectors[1,:]*self.reduced_Q[Q,1]+
                    lattice.r_lattice_vectors[2,:]*self.reduced_Q[Q,2])

    # ---------------------------------------------------------------------------------------------





