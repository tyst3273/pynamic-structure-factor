
# system modules
import numpy as np
import multiprocessing as mp
import itertools

# custom modules


class c_multi_processing:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,utils,comm):

        """
        tools for parallelism using multiprocessing
        """

        self.utils = utils
        self.comm = comm

        # number of processes to use for kpt parallelization
        self.num_kpt_procs = utils.config.num_kpt_procs
        if self.num_kpt_procs == 1:
            self.use_kpt_parallelism = False
            msg = 'parallism over k-points is disabled\n'
        else:
            self.use_kpt_parallelism = True    
            msg = f'number of processes for k-point parallelism: {self.num_kpt_procs}\n'

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

        _num_kpts = self.comm.k_points.num_kpts
        _num_procs = self.num_kpt_procs

        # the indices of the kpts on each process
        self.kpts_on_proc, self.num_kpts_per_proc = \
            self._distribute_round_robin(_num_kpts,_num_procs)

        msg_1 = 'number of k-points on each k-point process:'
        for ii in range(self.num_kpt_procs):

            _nk = self.num_kpts_per_proc[ii]

            _ = f'proc[{ii}]:'
            msg_1 = msg_1+f'\n  {_:9} {_nk}'

        # check if it will work
        if np.any(self.num_kpts_per_proc == 0):
            self.utils.log.crash('at least one process will treat 0 k-points!\n'\
                           ' decrease num_kpt_procs or use more k-points')

        # print wassup
        self.max_num_kpts_over_procs = self.num_kpts_per_proc.max()
        msg_2 = 'max number of k-points over all k-point processes:' \
                f' {self.max_num_kpts_over_procs}' 
        self.utils.log.info(msg_2,msg_label='k-point parallelism')

        if self.utils.config.verbosity > 0:
            self.utils.log.info(msg_1)

    # -----------------------------------------------------------------------------------------------

    def distribute_bonds_over_procs(self):

        """
        distribute the bonds over processes "round-robin" style
        """

        _num_hop = self.comm.hopping.num_bonds_to_calc
        _num_ovlp = self.comm.overlap.num_bonds_to_calc
        _num_procs = self.num_bond_procs

        # the indices of the hopping integrals on each process
        self.hop_on_proc, self.num_hop_per_proc = \
            self._distribute_round_robin(_num_hop,_num_procs)

        # the indices of the overlap integrals on each process
        self.ovlp_on_proc, self.num_ovlp_per_proc = \
            self._distribute_round_robin(_num_ovlp,_num_procs)

        msg_1 = 'number of integrals on each bond process:'
        for ii in range(_num_procs):

            _nh = self.num_hop_per_proc[ii]
            _no = self.num_ovlp_per_proc[ii]

            _ = f'proc[{ii}]:'
            msg_1 = msg_1+f'\n  {_:9}  hopping: {_nh}, overlap: {_no}'

        # check if it will work
        if np.any(self.num_hop_per_proc == 0) or np.any(self.num_ovlp_per_proc == 0):
            self.utils.log.crash('at least one process will treat 0 bonds!\n'\
                           ' decrease num_bond_procs')

        # print wassup
        self.max_num_bonds_over_procs = max(self.num_hop_per_proc.max(),
                        self.num_ovlp_per_proc.max())
        msg_2 = f'max number of bonds over all bond processes: {self.max_num_bonds_over_procs}'
        self.utils.log.info(msg_2,msg_label='bond parallelism')

        if self.utils.config.verbosity > 0:
            self.utils.log.info(msg_1)

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





