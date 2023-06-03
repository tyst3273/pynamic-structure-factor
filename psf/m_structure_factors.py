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
from scipy.fft import fft
import multiprocessing as mp

from psf.m_timing import _timer
from psf.m_error import crash


_thz2meV = 4.13566553853599

class c_structure_factors:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm,timers):

        """
        class to calculate structure factors from MD trajectories
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config
        self.timers = timers

        self.calc_sqw = config.calc_sqw

        _num_Q = self.comm.Qpoints.num_Q
        _num_steps = self.comm.traj.num_block_steps

        self.exp_iQr = np.zeros((_num_Q,_num_steps),dtype=complex)

        # always calculated
        self.sq_elastic = np.zeros(_num_Q,dtype=float)
        msg = '\n*** elastic intensity ***\ncalculating elastic intensity\n'
        msg += '  |<exp(iQ.r(t))>|**2 '
        print(msg)

        # optional
        if self.calc_sqw:
            self._setup_energies()
            self.sqw = np.zeros((_num_Q,self.num_energy),dtype=float)

            msg = '\n*** dynamic structure factor ***\n'
            msg += '  |FT[exp(iQ.r(t))]|**2 \n\n'
            msg += f'num_energy: {self.num_energy}\n'
            msg += f'energy_step:\n {self.energy_step: 8.6f} (meV)\n'
            msg += f' {self.energy_step/_thz2meV: 8.6f} (tHz)\n'
            msg += f'energy_max:\n {self.energy_max: 8.4f} (meV)\n'
            msg += f' {self.energy_max/_thz2meV: 8.4f} (tHz)'
            print(msg)

    # ----------------------------------------------------------------------------------------------

    def _setup_energies(self):

        """
        get energies based on size of blocks of traj
        """

        self.num_energy = self.comm.traj.num_block_steps

        # timestep is in fs, convert to freq. in tHz, then to meV
        self.energy = np.fft.fftfreq(self.num_energy, \
                self.comm.traj.effective_time_step*1e-3)*_thz2meV
        self.energy_step = self.energy[1]-self.energy[0]
        self.energy_max = self.energy.max()

    # ----------------------------------------------------------------------------------------------

    def calculate_structure_factors(self):

        """
        macro to calculate scattering
        """

        self.timers.start_timer('calc_strufacs',units='m')

        _num_blocks = self.comm.traj.num_block_avg
        msg = '\n*** structure factors ***\n'
        msg += f'calculating structure factors on {_num_blocks} blocks\n'
        if self.comm.paral.num_Qpoint_procs == 1:
            msg += 'doing calculation with 1 proc (serial)\n'
        else:
            msg += f'using {self.comm.paral.num_Qpoint_procs} processes (parallel)\n'
            msg += f'only printing info for proc-0\n'
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
        # and by number of steps to normalize vs traj. length
        _num_steps = self.comm.traj.num_block_steps
        _num_atoms = self.comm.traj.num_atoms

        self.sq_elastic /= _num_blocks*_num_atoms*_num_steps

        if self.calc_sqw:
            self.sqw /= _num_blocks*_num_atoms*_num_steps 

        self.timers.stop_timer('calc_strufacs')

    # ----------------------------------------------------------------------------------------------

    def _calculate_on_block(self,block):

        """
        reads data for a block, distributes data over procs, runs them, collects data, 
        writes it to appropriate arrays
        """

        self.timers.start_timer('calc_on_block',units='m')

        # get stuff for parallelization
        _num_procs = self.comm.paral.num_Qpoint_procs

        # get the trajectory from the file
        self.comm.traj.read_trajectory_block(block)

        # Queue for passing data between procs
        self.queue = mp.Queue()

        # do serial calculation if only 1 proc
        if _num_procs == 1:

            self._serial_loop_on_Q()

        # do parallel calculation 
        else:

            _procs = []
            for _proc in range(_num_procs):
                _procs.append(mp.Process(target=self._parallel_loop_on_Q,args=[_proc]))

            # start execution
            for _proc in _procs:
                _proc.start()

            # get the results from queue
            for _proc in range(_num_procs):

                proc, exp_iQr = self.queue.get()
                Q_inds = self.comm.paral.Qpts_on_proc[proc]
                self.exp_iQr[Q_inds,:] = exp_iQr[...]

            # close queue, kill it
            self.queue.close()
            self.queue.join_thread()

            # blocking; wait for all procs to finish before moving on
            for _proc in _procs:
                _proc.join()
                _proc.close()

        self._calc_strufacs(self.exp_iQr)

        self.timers.stop_timer('calc_on_block')

    # ----------------------------------------------------------------------------------------------

    def _calc_strufacs(self,exp_iQr):

        """
        take exp_iQr (rho(Q,t)) and calculate everything from it
        """

        self.sq_elastic += np.abs(np.mean(exp_iQr,axis=1))**2

        if self.calc_sqw:
            self.sqw += np.abs(fft(exp_iQr,axis=1,norm='forward'))**2

    # ----------------------------------------------------------------------------------------------

    def _parallel_loop_on_Q(self,proc=0):

        """
        loops over the Q-points assinged to the proc, calcs rho(Q,t), then puts data in queue
        """

        # convenience refs for below
        _exp_type = self.comm.xlengths.experiment_type
        _num_steps = self.comm.traj.num_block_steps
        _num_atoms = self.comm.traj.num_atoms

        # get the Q-points that this proc is supposed do
        _Qpt_inds = self.comm.paral.Qpts_on_proc[proc]
        _Qpts = self.comm.Qpoints.Q_cart[_Qpt_inds,:]
        _nQ = _Qpts.shape[0]

        exp_iQr = np.zeros((_nQ,_num_steps),dtype=complex)

        # get this for neutrons once and for all
        if _exp_type == 'neutrons':
            _x = self.comm.xlengths.scattering_lengths

        # only print info for proc-0
        if proc == 0:
            msg = f'there are {_nQ} Q-points to do ...'
            print(msg)

        # now loop over Q-points on this proc
        for ii in range(_nQ):

            if proc == 0:
                msg = f'  now on Qpt[{ii}]'
                print(msg)

            # depends on Q for xrays, but calculated earlier (only looked up here)
            if _exp_type == 'xrays':
                _Q_ind = _Qpt_inds[ii]
                _x = self.comm.xlengths.form_factors[:,_Q_ind]

            # vectorized Q.r dot product, sum over atoms gives space FT
            _Q = _Qpts[ii,:].reshape(1,3)
            _x_tile = np.tile(_x.reshape(1,_num_atoms),reps=(_num_steps,1))
            _exp_iQr = np.tile(_Q,reps=(_num_steps,_num_atoms,1))
            _exp_iQr = np.sum(_exp_iQr*self.comm.traj.pos,axis=2)
            _exp_iQr = np.exp(1j*_exp_iQr)*_x_tile
            _exp_iQr = np.sum(_exp_iQr,axis=1)

            exp_iQr[ii,:] = _exp_iQr

        # put results in queue to be passed to main proc
        self.queue.put([proc,exp_iQr])
        
    # ---------------------------------------------------------------------------------------------- 
    
    def _serial_loop_on_Q(self):

        """
        loops over all of the Q-points on a single proc, calcs rho(Q,t) 
        """

        # convenience refs for below
        _exp_type = self.comm.xlengths.experiment_type
        _num_steps = self.comm.traj.num_block_steps
        _num_atoms = self.comm.traj.num_atoms

        # get the Q-points that this proc is supposed do
        _Qpt_inds = self.comm.paral.Qpts_on_proc[0]
        _Qpts = self.comm.Qpoints.Q_cart[_Qpt_inds,:]
        _nQ = _Qpts.shape[0]

        #exp_iQr = np.zeros((_nQ,_num_steps),dtype=complex)

        # get this for neutrons once and for all
        if _exp_type == 'neutrons':
            _x = self.comm.xlengths.scattering_lengths

        msg = f'there are {_nQ} Q-points to do ...'
        print(msg)

        # now loop over Q-points on this proc
        for ii in range(_nQ):

            msg = f'  now on Qpt[{ii}]'
            print(msg)

            _Q_ind = _Qpt_inds[ii]

            # depends on Q for xrays, but calculated earlier (only looked up here)
            if _exp_type == 'xrays':
                _x = self.comm.xlengths.form_factors[:,_Q_ind]

            # vectorized Q.r dot product, sum over atoms gives space FT
            _Q = _Qpts[ii,:].reshape(1,3)
            _x_tile = np.tile(_x.reshape(1,_num_atoms),reps=(_num_steps,1))
            _exp_iQr = np.tile(_Q,reps=(_num_steps,_num_atoms,1))
            _exp_iQr = np.sum(_exp_iQr*self.comm.traj.pos,axis=2)
            _exp_iQr = np.exp(1j*_exp_iQr)*_x_tile
            _exp_iQr = np.sum(_exp_iQr,axis=1)

            self.exp_iQr[_Q_ind,:] = _exp_iQr

    # ---------------------------------------------------------------------------------------------- 

    def put_on_mesh(self):

        """
        convenience methods that put S(Q,w) on list of Q onto cartesian grid of Q-points

        this calls methods from c_Qpoints
        """

        _use_mesh = self.comm.Qpoints.use_mesh
        if not _use_mesh:
            msg = 'a Qpoint-mesh was not requested, i.e. cannot be unfolded!\n' \
                  'use \'Qpoints_option\' = \'mesh\' or \'write_mesh\' \n'
            crash(msg)

        if self.calc_sqw:
            self.sqw = self.comm.Qpoints.unfold_onto_Q_mesh(self.sqw)

        self.sq_elastic = self.comm.Qpoints.unfold_onto_Q_mesh(self.sq_elastic)

    # ----------------------------------------------------------------------------------------------






