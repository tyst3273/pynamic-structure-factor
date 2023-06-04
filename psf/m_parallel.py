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

# system modules
import numpy as np
import multiprocessing as mp
import itertools

# custom modules
from psf.m_error import crash

class c_multi_processing:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm):

        """
        tools for parallelism using multiprocessing
        """

        self.config = config
        self.comm = comm

        # number of processes to use for kpt parallelization
        self.num_Qpoint_procs = config.num_Qpoint_procs

        # number of processes to use for bond calc. parallelization
        self.distribute_Q_over_procs()

    # ----------------------------------------------------------------------------------------------

    def distribute_Q_over_procs(self):

        """
        distribute the kpts over processes "round-robin" style
        """

        num_Qpts = self.comm.Qpoints.num_Q
        num_procs = self.num_Qpoint_procs

        # the indices of the kpts on each process
        self.Qpts_on_proc, self.num_Qpts_per_proc = \
            self._distribute_round_robin(num_Qpts,num_procs)

        # check if it will work
        if np.any(self.num_Qpts_per_proc == 0):
            crash('at least one process will treat 0 Q-points!\n'\
                  'decrease num_Qpoint_procs or use more Q-points\n')

        # print wassup
        self.max_num_Qpts_on_procs = self.num_Qpts_per_proc.max()
        msg = 'max number of Q-points over all processes:' \
                f' {self.max_num_Qpts_on_procs}\n'

        msg += 'number of Q-points on each Q-point process:'
        for ii in range(self.num_Qpoint_procs):

            _nQ = self.num_Qpts_per_proc[ii]
            _ = f'proc[{ii}]:'
            msg = msg+f'\n  {_:9} {_nQ}'

        print('\n*** Q-point parallelism ***')
        print(msg)

    # -----------------------------------------------------------------------------------------------

    def _distribute_round_robin(self,num_inds,num_procs):

        """
        distribute inds over procs round-robin style (done automatically with numpy)
        """

        _on_proc = np.array_split(np.arange(num_inds),num_procs)
        _num_per_proc = [len(_) for _ in _on_proc]
        _num_per_proc = np.array(_num_per_proc)

        return _on_proc, _num_per_proc

    # -----------------------------------------------------------------------------------------------


