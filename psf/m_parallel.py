
# system modules
import numpy as np
import multiprocessing as mp
import itertools

# custom modules


class c_multi_processing:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,comm):

        """
        tools for parallelism using multiprocessing
        """

        # number of processes to use for kpt parallelization
        self.num_Qpoint_procs = utils.config.num_Qpoint_procs
        if self.num_Qpoint_procs == 1:
            self.use_kpt_parallelism = False
            msg = 'parallism over k-points is disabled\n'
        else:
            self.use_kpt_parallelism = True    
            msg = f'number of processes for k-point parallelism: {self.num_Qpoint_procs}\n'

        # number of processes to use for bond calc. parallelization
        self.num_bond_procs = utils.config.num_bond_procs
        if self.num_bond_procs == 1:
            self.use_bond_parallelism = False
            msg = msg+' parallism over bond calculations disabled'
        else:
            self.use_bond_parallelism = True
            msg = msg+f' number of processes for bond parallelism: {self.num_bond_procs}'

        self.utils.log.console(msg)
        self.utils.log.info(msg,msg_label='parallelism')

    # ----------------------------------------------------------------------------------------------

    def distribute_k_over_procs(self):

        """
        distribute the kpts over processes "round-robin" style
        """

        _num_Qpts = self.comm.k_points.num_kpts
        _num_procs = self.num_Qpoint_procs

        # the indices of the kpts on each process
        self.kpts_on_proc, self.num_Qpts_per_proc = \
            self._distribute_round_robin(_num_kpts,_num_procs)

        msg = 'number of Q-points on each Q-point process:'
        for ii in range(self.num_Qpoint_procs):

            _nk = self.num_Qpts_per_proc[ii]

            _ = f'proc[{ii}]:'
            msg = msg+f'\n  {_:9} {_nk}'

        # check if it will work
        if np.any(self.num_Qpts_per_proc == 0):
            self.utils.log.crash('at least one process will treat 0 Q-points!\n'\
                           ' decrease num_Qpoint_procs or use more Q-points')

        # print wassup
        self.max_num_Qpts_on_procs = self.num_Qpts_per_proc.max()
        msg = 'max number of Q-points over all processes:' \
                f' {self.max_num_Qpts_on_procs}' 

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





