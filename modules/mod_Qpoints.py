#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                                           !
#   ! Copyright 2021 by Tyler C. Sterling and Dmitry Reznik,                    !
#   ! University of Colorado Boulder                                            !
#   !                                                                           !
#   ! This file is part of the pynamic-structure-factor (PSF) software.         !
#   ! PSF is free software: you can redistribute it and/or modify it under      !
#   ! the terms of the GNU General Public License as published by the           !
#   ! Free software Foundation, either version 3 of the License, or             !
#   ! (at your option) any later version. PSF is distributed in the hope        !
#   ! that it will be useful, but WITHOUT ANY WARRANTY; without even the        !
#   ! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  !
#   ! See the GNU General Public License for more details.                      !
#   !                                                                           !
#   ! A copy of the GNU General Public License should be available              !
#   ! alongside this source in a file named gpl-3.0.txt. If not see             !
#   ! <http://www.gnu.org/licenses/>.                                           !
#   !                                                                           !
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import numpy as np
from modules.mod_utils import print_stdout, PSF_exception

class Qpoints:

    """
    store, intialize, and distribute Qpoints over procs
    """

    # -----------------------------------------------------------------------------------------

    def generate_Qpoints(self,invars,lattice):

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

            if Q_count == 49: # break if >= 50 Q points
                message = ('...............................\n number of Qpoints is >= 50.'
                           '\n suppressing output.')
                print_stdout(message)
                break
            Q_count = Q_count+1

        self._convert_Q_to_1_over_Angstrom(lattice)      # convert the Q from rlu to 1/A

    # --------------------------------------------------------------------------------------------

    def distribute_Q_over_procs(self,num_procs):

        """
        num_Q_per_proc is determined as the largest integer dividing the total_Qsteps number. 
        the remainder is placed on rank 0 (if there is a remainder...)
        """

        # set up array of how many Q on each proc
        nQpp_arr = np.zeros(num_procs).astype(int) # Num_Q_Per_Processor_ARRay
        proc = 0
        for qq in range(self.total_Qsteps):
            if proc == num_procs:
                proc = 0
            nQpp_arr[proc] = nQpp_arr[proc]+1
            proc = proc+1

        # if any procs have 0 Q, print error message and exit
        if nQpp_arr.min() == 0:
            message = 'atleast one processor will do 0 Qpoints.\n' \
                           ' increase number of procs or decrease number of Q points'
            raise PSF_exception(message)

        # print parellelism info
        message = f'process: 0    Q points: {nQpp_arr[0]:g}\n'
        for ii in range(1,num_procs):
            message = message+f' process: {ii:g}    Q points: {nQpp_arr[ii]:g}\n'
        print_stdout(message,msg_type='Q points on each process')

        # put the Qpoint indicies for each proc in the list
        self.Q_on_procs = []
        shift = 0
        for ii in range(num_procs):
            self.Q_on_procs.append(list(range(shift,shift+nQpp_arr[ii])))
            shift = shift+nQpp_arr[ii]

    # ---------------------------------------------------------------------------------------------

    def reconvert_Q_points(self,lattice):

        """
        interface to use if recalcuting cell lengths from traj file: the recip. lattice changes
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

        self.total_Qpoints = np.zeros((self.total_Qsteps,3))
        for Q in range(self.total_Qsteps):
            self.total_Qpoints[Q,:] = (lattice.r_lattice_vectors[0,:]*self.total_reduced_Q[Q,0]+ 
                    lattice.r_lattice_vectors[1,:]*self.total_reduced_Q[Q,1]+ 
                    lattice.r_lattice_vectors[2,:]*self.total_reduced_Q[Q,2])

    # ---------------------------------------------------------------------------------------------





