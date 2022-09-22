
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

        _num_Qpts = self.comm.Qpoints.num_Q
        _num_procs = self.num_Qpoint_procs

        # the indices of the kpts on each process
        self.Qpts_on_proc, self.num_Qpts_per_proc = \
            self._distribute_round_robin(_num_Qpts,_num_procs)

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
        distribute inds over procs round-robin style
        """

        # the indices on each process
        _on_proc = []
        for ii in range(num_procs):
            _on_proc.append([])

        # iterator for the procesors
        _iter = itertools.cycle(range(num_procs))

        # assign inds to each process round robin style
        for ii in range(num_inds):
            _on_proc[(next(_iter))].append(ii)

        # number of kpts for each process to do
        _num_per_proc = np.zeros(num_procs,dtype=int)

        # count number of inds on each proc
        for ii in range(num_procs):

            _num = len(_on_proc[ii])
            _num_per_proc[ii] = _num

        return _on_proc, _num_per_proc

    # -----------------------------------------------------------------------------------------------





