
import numpy as np
from scipy.fft import fft
from psf.m_timing import _timer


_thz2meV = 4.13566553853599

class c_structure_factors:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm,timers):

        """
        template
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config
        self.timers = timers

        self.calc_sqw = config.calc_sqw
        self.calc_diffuse = config.calc_diffuse
        self.calc_bragg = config.calc_bragg

        # set stuff up
        if self.calc_sqw:
            self._setup_energies()
            self.sqw = np.zeros((self.comm.Qpoints.num_Q,self.num_energy),dtype=float)

            msg = '\n*** dynamic structure factor ***\ncalculating S(Q,E) ~ |rho(Q,E)|**2 \n'
            msg += f'num_energy: {self.num_energy}\n'
            msg += f'energy_step:\n {self.energy_step: 8.6f} (meV)\n'
            msg += f' {self.energy_step/_thz2meV: 8.6f} (tHz)\n'
            msg += f'energy_max:\n {self.energy_max: 8.4f} (meV)\n'
            msg += f' {self.energy_max/_thz2meV: 8.4f} (tHz)'
            print(msg)

        if self.calc_diffuse:   
            self.fqt_sq = np.zeros((self.comm.Qpoints.num_Q,self.comm.traj.num_block_steps),dtype=float)
            self.sq_diffuse = np.zeros(self.comm.Qpoints.num_Q,dtype=float)

            msg = '\n*** diffuse intensity ***\ncalculating diffuse |F(Q,t)|**2 and \n'
            msg += '  S(Q) ~ <|F(Q,t)|**2> '
            print(msg)

        if self.calc_bragg:
            self.sq_bragg = np.zeros(self.comm.Qpoints.num_Q,dtype=float)

            msg = '\n*** bragg intensity ***\ncalculating bragg intensity\n'
            msg += '  I(Q) ~ |<F(Q),t)>|**2 '
            print(msg)

    # ----------------------------------------------------------------------------------------------

    def _setup_energies(self):

        """
        get energies based on size of blocks of traj
        """

        self.num_energy = self.comm.traj.num_block_steps

        # timestep is in fs, convert to freq. in tHz, then to meV
        self.energy = np.fft.fftfreq(self.num_energy,self.comm.traj.md_time_step*1e-3)*_thz2meV
        self.energy_step = self.energy[1]-self.energy[0]
        self.energy_max = self.energy.max()

    # ----------------------------------------------------------------------------------------------

    def calculate_structure_factors(self):

        """
        macro within c_structure_factors to calculate scattering intensity
        """

        self.timers.start_timer('calc_strufacs',units='m')

        _num_blocks = self.comm.traj.num_block_avg
        msg = '\n*** structure factors ***\n'
        msg += f'calculating structure factors on {_num_blocks} blocks '
        msg += f'using {self.comm.paral.num_Qpoint_procs} processes\n'
        msg += f'only printing info for proc-0\n\n'
        msg += 'proc-0 treats >= the number of Q-points on other processes\n'
        msg += f'number of Q-points on proc-0: {self.comm.paral.max_num_Qpts_on_procs}\n'
        print(msg)

        for ii, _block in enumerate(self.comm.traj.blocks):
        
            _ = f'block[{ii}]:'

            # start a local timer
            _block_timer = _timer(_,units='s')
            _block_timer.start()

            msg = '-------------------------------------------------------------\n\n'
            msg += f'{_:12} {_block}'
            print(msg)

            # calculate for this block; internally parallelized
            self._calculate_on_block(_block)

            # end the local timer
            _block_timer.stop()
            _block_timer.print_timing()

        # divide by number of blocks and number of atoms (to normalize vs system size)
        _num_atoms = self.comm.traj.num_atoms
        if self.calc_bragg:
            self.sq_bragg /= _num_blocks #*_num_atoms
        if self.calc_diffuse:
            self.sq_bragg /= _num_blocks #*_num_atoms
        if self.calc_sqw:
            self.sq_bragg /= _num_blocks #*_num_atoms

        self.timers.stop_timer('calc_strufacs')

    # ----------------------------------------------------------------------------------------------

    def _calculate_on_block(self,_block):

        """
        reads data for a block, distributes data over procs, runs them, collects data, 
        writes it to appropriate arrays
        """

        self.timers.start_timer('calc_on_block',units='m')

        # get the trajectory from the file
        self.comm.traj.read_trajectory_block(_block)

        # place holder until i set up multiprocessing here!

        if True:

            # get the results
            res = self._proc_loop_on_Q(proc=0) 

            # get the proc
            proc = res[0]
            res.pop(0)

            # get the Q-point indices this proc did
            Q_inds = self.comm.paral.Qpts_on_proc[proc]

            # get the data and put on correct Qpts
            if self.calc_bragg:
                self.sq_bragg[Q_inds] += res[0]
                res.pop(0)
            if self.calc_diffuse:
                self.fqt_sq[Q_inds,:] += res[0]
                res.pop(0)
                self.sq_diffuse[Q_inds] += res[0]
                res.pop(0)
            if self.calc_sqw:
                self.sqw[Q_inds,:] += res[0]
                res.pop(0)

        self.timers.stop_timer('calc_on_block')

    # ----------------------------------------------------------------------------------------------

    def _proc_loop_on_Q(self,proc=0):

        """
        loops over the Q-points assinged to the proc, calcs S(Q,E) etc, the puts data in queue
        """

        self.timers.start_timer('loop_on_Q',units='s')

        # convenience refs for below
        _exp_type = self.comm.xlengths.experiment_type
        _num_steps = self.comm.traj.num_block_steps
        _num_atoms = self.comm.traj.num_atoms
        _num_energy = self.num_energy

        # get the Q-points that this proc is supposed do
        _Qpt_inds = self.comm.paral.Qpts_on_proc[proc]
        _Qpts = self.comm.Qpoints.Q_cart[_Qpt_inds,:]
        _nQ = _Qpts.shape[0]

        # set stuff up
        if self.calc_sqw:
            _sqw = np.zeros((_nQ,_num_energy),dtype=float)
        if self.calc_diffuse:
            _fqt_sq = np.zeros((_nQ,_num_steps),dtype=float)
            _sq_diffuse = np.zeros(_nQ,dtype=float)
        if self.calc_bragg:
            _sq_bragg = np.zeros(_nQ,dtype=float)

        # get this for neutrons once and for all
        if _exp_type == 'neutrons':
            _x = self.comm.xlengths.scattering_lengths

        if proc == 0:
            msg = f'there are {_nQ} Q-points to do ...'
            print(msg)

        # now loop over Q-points on this proc
        for ii in range(_nQ):

            if proc == 0:
                msg = f'  now on Qpt[{ii}]'
                print(msg)

            _Q = _Qpts[ii,:].reshape(1,3)
            _Q_ind = _Qpt_inds[ii]

            # depends on Q for xrays, but calculate earlier
            if _exp_type == 'xrays':
                _x = self.comm.xlengths.form_factors[:,_Q_ind]

            # vectorized Q.r dot product, sum over atoms gives space FT
            _x = np.tile(_x.reshape(1,_num_atoms),reps=(_num_steps,1))
            _exp_iQr = np.tile(_Q,reps=(_num_steps,_num_atoms,1))
            _exp_iQr = np.sum(_exp_iQr*self.comm.traj.pos,axis=2)
            _exp_iQr = np.exp(1j*_exp_iQr)*_x
            _exp_iQr = np.sum(_exp_iQr,axis=1)

            # calculate whatever was requested
            if self.calc_bragg:
                _sq_bragg[ii] = np.abs((_exp_iQr).mean())**2
            if self.calc_diffuse:
                _fqt_sq[ii,:] = np.abs(_exp_iQr)**2
                _sq_diffuse[ii] = _fqt_sq[ii,:].mean()
            if self.calc_sqw:
                _sqw[ii,:] = np.abs(fft(_exp_iQr))**2

        # return list of results
        res = [proc]
        if self.calc_bragg:
            res.append(_sq_bragg)
        if self.calc_diffuse:
            res.append(_fqt_sq)
            res.append(_sq_diffuse)
        if self.calc_sqw:
            res.append(_sqw)

        self.timers.stop_timer('loop_on_Q')
        
        return res

    # ---------------------------------------------------------------------------------------------- 

