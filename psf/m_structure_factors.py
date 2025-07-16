#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                                           !
#   ! Copyright 2025 by Tyler C. Sterling                                       !
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
import multiprocess as mp

from psf.m_timing import _timer
from psf.m_error import crash


_switch = 10 # how often to print Q timer
_thz2meV = 4.13566553853599

# --------------------------------------------------------------------------------------------------

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
        self.calc_coherent = config.calc_coherent
        self.calc_incoherent = config.calc_incoherent

        _num_Q = self.comm.Qpoints.num_Q
        _num_steps = self.comm.traj.num_block_steps

        self.exp_iQr = np.zeros((_num_Q,_num_steps),dtype=complex)

        # elastic part always calculated
        msg = '\n*** elastic intensity ***'
        if self.calc_coherent:
            self.coh_sq_elastic = np.zeros(_num_Q,dtype=float)
            msg += '\ncalculating coherent elastic intensity'    
        if self.calc_incoherent:
            self.inc_sq_elastic = np.zeros(_num_Q,dtype=float)
            msg += '\ncalculating incoherent elastic intensity'
        print(msg)

        # optional
        if self.calc_sqw:

            self._setup_energies()

            msg = '\n*** dynamic structure factor ***\n'
            msg += f'num_energy: {self.num_energy}\n'
            msg += f'energy_step:\n {self.energy_step: 8.6f} (meV)\n'
            msg += f' {self.energy_step/_thz2meV: 8.6f} (tHz)\n'
            msg += f'energy_max:\n {self.energy_max: 8.4f} (meV)\n'
            msg += f' {self.energy_max/_thz2meV: 8.4f} (tHz)\n'

            if self.calc_coherent:
                self.coh_sqw = np.zeros((_num_Q,self.num_energy),dtype=float)
                msg += '\ncalculating coherent inelastic intensity'
            if self.calc_incoherent:
                self.inc_sqw = np.zeros((_num_Q,self.num_energy),dtype=float)
                msg += '\ncalculating incoherent inelastic intensity'
            
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
        if self.comm.paral.num_Qpoint_procs > 1:
            msg += f'using {self.comm.paral.num_Qpoint_procs} processes (parallel)\n'
            msg += f'only printing info for proc-0\n'
            msg += 'proc-0 treats >= the number of Q-points on other processes\n'
            msg += f'number of Q-points on proc-0: {self.comm.paral.max_num_Qpts_on_procs}\n'
        else:
            msg += f'using 1 processe (serial)\n'
            msg += f'number of Q-points: {self.comm.Qpoints.num_Q}\n'
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

        msg = '-------------------------------------------------------------'
        print(msg)

        _num_steps = self.comm.traj.num_block_steps
        _num_atoms = self.comm.traj.num_atoms

        # average over blocks
        if self.calc_coherent:
            self.coh_sq_elastic /= _num_blocks * _num_steps * _num_atoms
        if self.calc_incoherent:
            self.inc_sq_elastic /= _num_blocks * _num_steps * _num_atoms

        if self.calc_sqw:
            if self.calc_coherent:
                self.coh_sqw /= _num_blocks * _num_steps * _num_atoms
            if self.calc_incoherent:
                self.inc_sqw /= _num_blocks * _num_steps * _num_atoms

        self.timers.stop_timer('calc_strufacs')

    # ----------------------------------------------------------------------------------------------

    def _calculate_on_block(self,block):

        """
        reads data for a block, distributes data over procs, runs them, collects data, 
        writes it to appropriate arrays
        """

        self.timers.start_timer('calc_on_block',units='m')

        _num_procs = self.comm.paral.num_Qpoint_procs

        # get the trajectory from the file
        self.comm.traj.read_trajectory_block(block)

        # do serial calculation if only 1 proc
        if _num_procs == 1:

            self._serial_loop_on_Q()

        else:
        
            # Queue for passing data between procs
            self.queue = mp.Queue()

            _procs = []
            for _proc in range(_num_procs):
                _procs.append(mp.Process(target=self._parallel_loop_on_Q,args=[_proc]))

            # start execution
            for _proc in _procs:
                _proc.start()

            # get the results from queue
            for _proc in range(_num_procs):

                self.timers.start_timer('get_from_queue',units='s')

                _queue = self.queue.get()

                proc = _queue.pop(0)
                Q_inds = self.comm.paral.Qpts_on_proc[proc]

                # same order as put into queue
                if self.calc_coherent:
                    _coh_sq_elastic = _queue.pop(0)
                    self.coh_sq_elastic[Q_inds] += _coh_sq_elastic  
                if self.calc_incoherent:
                    _inc_sq_elastic = _queue.pop(0)
                    self.inc_sq_elastic[Q_inds] += _inc_sq_elastic
                if self.calc_sqw:
                    if self.calc_coherent:
                        _coh_sqw = _queue.pop(0)
                        self.coh_sqw[Q_inds,:] += _coh_sqw
                    if self.calc_incoherent:
                        _inc_sqw = _queue.pop(0)
                        self.inc_sqw[Q_inds,:] += _inc_sqw
                
                if proc == 0:
                    _timers = _queue.pop(0)
                    for _key in _timers.timers:
                        if _key not in self.timers.timers:
                            self.timers.timers[_key] = _timers.timers[_key]

                self.timers.stop_timer('get_from_queue')

            # close queue, kill it
            self.queue.close()
            self.queue.join_thread()

            # blocking; wait for all procs to finish before moving on
            for _proc in _procs:
                _proc.join()
                _proc.close()

        self.timers.stop_timer('calc_on_block')

    # ----------------------------------------------------------------------------------------------

    def _parallel_loop_on_Q(self,proc=0):

        """
        loops over the Q-points assinged to the proc, calcs rho(Q,t), then puts data in queue
        """

        # convenience refs for below
        _exp_type = self.comm.xlengths.experiment_type
        _num_steps = self.comm.traj.num_block_steps
        _num_atoms = self.comm.traj.num_atoms

        _num_fft_procs = self.config.num_fft_procs

        # get the Q-points that this proc is supposed do
        _Qpt_inds = self.comm.paral.Qpts_on_proc[proc]
        _Qpts = self.comm.Qpoints.Q_cart[_Qpt_inds,:]
        _nQ = _Qpts.shape[0]

        if self.calc_coherent:
            coh_sq_elastic = np.zeros(_nQ,dtype=float)
        if self.calc_incoherent:
            inc_sq_elastic = np.zeros(_nQ,dtype=float)

        if self.calc_sqw:
            _num_energy = self.num_energy
            if self.calc_coherent:
                coh_sqw = np.zeros((_nQ,_num_energy),dtype=float)
            if self.calc_incoherent:
                inc_sqw = np.zeros((_nQ,_num_energy),dtype=float)

        # get this for neutrons once and for all
        if _exp_type == 'neutrons':
            if self.calc_coherent:
                _b = self.comm.xlengths.coherent_scattering_lengths
                _b_tile = np.tile(_b.reshape(1,_num_atoms),reps=(_num_steps,1)) 
            if self.calc_incoherent:
                _x = self.comm.xlengths.incoherent_scattering_cross_section
                _x_tile = np.tile(_x.reshape(_num_atoms,1),reps=(1,_num_steps)) 

        # only print info for proc-0
        if proc == 0:
            msg = f'there are {_nQ} Q-points to do ...'
            print(msg)

        # stuff for intermittent timing
        if proc == 0:
            _Q_timer = _timer('Q_timer',units='s')

        # now loop over Q-points on this proc
        for ii in range(_nQ):

            if proc == 0:
                
                msg = f'  now on Qpt[{ii}]'
                print(msg)

                # start a local timer
                _Q_timer.start()

            # depends on Q for xrays, but calculated earlier. only coherent scattering for xrays
            if _exp_type == 'xrays':
                _Q_ind = _Qpt_inds[ii]
                _b = self.comm.xlengths.form_factors[:,_Q_ind]
                _b_tile = np.tile(_b.reshape(1,_num_atoms),reps=(_num_steps,1)) 

            self.timers.start_timer('reshape',units='s')
            _Q = _Qpts[ii,:].reshape(1,3)
            _exp_iQr = np.tile(_Q,reps=(_num_steps,_num_atoms,1)) 
            _exp_iQr = np.sum(_exp_iQr*self.comm.traj.pos,axis=2) # Q.r
            _exp_iQr = np.exp(1j*_exp_iQr)
            self.timers.stop_timer('reshape')

            # calculate coherent scattering
            if self.calc_coherent:

                self.timers.start_timer('calc_coherent',units='s')

                # sum over atoms => space FT
                _F_Qt = np.sum(_b_tile*_exp_iQr,axis=1)            

                # integrate over time
                coh_sq_elastic[ii] = np.abs(np.mean(_F_Qt))**2

                # do time FT
                if self.calc_sqw:

                    self.timers.start_timer('coherent_FFT',units='s')
                    coh_sqw[ii,:] = np.abs(fft(_F_Qt,norm='forward'))**2
                    self.timers.stop_timer('coherent_FFT')
                
                self.timers.stop_timer('calc_coherent')

            # calculate incoherent scattering
            if self.calc_incoherent:

                self.timers.start_timer('calc_incoherent',units='s')

                # we want time in the rows (C order)
                _exp_iQr = _exp_iQr.T

                # integrate over time and sum over atoms
                inc_sq_elastic[ii] = np.sum(_x * np.abs(np.mean(_exp_iQr,axis=1))**2 )

                # do time FT
                if self.calc_sqw:

                    self.timers.start_timer('incoherent_FFT',units='s')
                    inc_sqw[ii,:] = np.sum( _x_tile * np.abs(fft(_exp_iQr,norm='forward',axis=1,
                                                            workers=_num_fft_procs))**2, axis=0 )
                    self.timers.stop_timer('incoherent_FFT')
                    
                self.timers.stop_timer('calc_incoherent')

            if proc == 0:
                
                _Q_timer.stop()

                if ii % _switch == 0 and ii != 0:    
                    _Q_timer.print_timing()

        # put results in queue to be passed to main proc
        _queue = [proc]
        if self.calc_coherent:
            _queue.append(coh_sq_elastic)
        if self.calc_incoherent:
            _queue.append(inc_sq_elastic)
        if self.calc_sqw:
            if self.calc_coherent:
                _queue.append(coh_sqw)
            if self.calc_incoherent:
                _queue.append(inc_sqw)

        if proc == 0:
            _queue.append(self.timers)

        self.queue.put(_queue)
    
    # ----------------------------------------------------------------------------------------------

    def _serial_loop_on_Q(self):

        """
        loops over the Q-points assinged to the proc, calcs rho(Q,t), then puts data in queue
        """

        # convenience refs for below
        _exp_type = self.comm.xlengths.experiment_type
        _num_steps = self.comm.traj.num_block_steps
        _num_atoms = self.comm.traj.num_atoms

        _num_fft_procs = self.config.num_fft_procs

        # get the Q-points that this proc is supposed do
        _Qpts = self.comm.Qpoints.Q_cart
        _nQ = _Qpts.shape[0]

        # get this for neutrons once and for all
        if _exp_type == 'neutrons':
            if self.calc_coherent:
                _b = self.comm.xlengths.coherent_scattering_lengths
                _b_tile = np.tile(_b.reshape(1,_num_atoms),reps=(_num_steps,1)) 
            if self.calc_incoherent:
                _x = self.comm.xlengths.incoherent_scattering_cross_section
                _x_tile = np.tile(_x.reshape(_num_atoms,1),reps=(1,_num_steps)) 

        msg = f'there are {_nQ} Q-points to do ...'
        print(msg)

        # stuff for intermittent timing
        _Q_timer = _timer('Q_timer',units='s')

        # now loop over Q-points on this proc
        for ii in range(_nQ):
            
            msg = f'  now on Qpt[{ii}]'
            print(msg)

            # start a local timer
            _Q_timer.start()

            # depends on Q for xrays, but calculated earlier. only coherent scattering for xrays
            if _exp_type == 'xrays':
                _b = self.comm.xlengths.form_factors[:,ii]
                _b_tile = np.tile(_b.reshape(1,_num_atoms),reps=(_num_steps,1)) 

            self.timers.start_timer('reshape',units='s')

            _Q = _Qpts[ii,:].reshape(1,3)
            _exp_iQr = np.tile(_Q,reps=(_num_steps,_num_atoms,1)) 
            _exp_iQr = np.sum(_exp_iQr*self.comm.traj.pos,axis=2) # Q.r
            _exp_iQr = np.exp(1j*_exp_iQr)

            self.timers.stop_timer('reshape')

            # calculate coherent scattering
            if self.calc_coherent:

                self.timers.start_timer('calc_coherent',units='s')

                # sum over atoms => space FT
                _F_Qt = np.sum(_b_tile*_exp_iQr,axis=1)            

                # integrate over time
                self.coh_sq_elastic[ii] += np.abs(np.mean(_F_Qt))**2

                # do time FT
                if self.calc_sqw:

                    self.timers.start_timer('coherent_FFT',units='s')
                    self.coh_sqw[ii,:] += np.abs(fft(_F_Qt,norm='forward'))**2
                    self.timers.stop_timer('coherent_FFT')

                self.timers.stop_timer('calc_coherent')

            # calculate incoherent scattering
            if self.calc_incoherent:

                self.timers.start_timer('calc_incoherent',units='s')

                # we want time in the rows (C order)
                _exp_iQr = _exp_iQr.T

                # integrate over time and sum over atoms
                self.inc_sq_elastic[ii] += np.sum(_x * np.abs(np.mean(_exp_iQr,axis=1))**2 )

                # do time FT
                if self.calc_sqw:

                    self.timers.start_timer('incoherent_FFT',units='s')
                    self.inc_sqw[ii,:] += np.sum( _x_tile * np.abs(fft(_exp_iQr,norm='forward',
                                                    axis=1,workers=_num_fft_procs))**2, axis=0 )
                    self.timers.stop_timer('incoherent_FFT')

                self.timers.stop_timer('calc_incoherent')
            
            _Q_timer.stop()

            if ii % _switch == 0 and ii != 0:    
                _Q_timer.print_timing()

    # ---------------------------------------------------------------------------------------------- 

    def unfold_structure_factors(self):

        """
        if reduced onto irreducible set, unfold onto full set. moreover, if on a mesh, unfold from 
        [num_Q]x[3] to [num_H]x[num_K]x[num_L] mesh
        """

        _calc_sqw = self.calc_sqw
        _calc_coh = self.calc_coherent
        _calc_inc = self.calc_incoherent
        _use_mesh = self.comm.Qpoints.use_mesh
        _Qpts = self.comm.Qpoints
        _symmetry = self.comm.symmetry

        # if using symmetry, put back onto full set
        if _symmetry.use_symmetry:

            if _calc_coh:
                self.coh_sq_elastic = _symmetry.unfold_onto_full_Q_set(self.coh_sq_elastic)
            if _calc_inc:
                self.inc_sq_elastic = _symmetry.unfold_onto_full_Q_set(self.inc_sq_elastic)

            if _calc_sqw:
                if _calc_coh:
                    self.coh_sqw = _symmetry.unfold_onto_full_Q_set(self.coh_sqw)
                if _calc_inc:
                    self.inc_sqw = _symmetry.unfold_onto_full_Q_set(self.inc_sqw)

        # if on Cartesian mesh, unfold from [num_Q]x[3] array to [num_H]x[num_K]x[num_L] mesh
        if _use_mesh:

            if _calc_coh:
                self.coh_sq_elastic = _Qpts.unfold_onto_cartesian_mesh(self.coh_sq_elastic)
            if _calc_inc:
                self.inc_sq_elastic = _Qpts.unfold_onto_cartesian_mesh(self.inc_sq_elastic)

            if _calc_sqw:
                if _calc_coh:
                    self.coh_sqw = _Qpts.unfold_onto_cartesian_mesh(self.coh_sqw)
                if _calc_inc:
                    self.inc_sqw = _Qpts.unfold_onto_cartesian_mesh(self.inc_sqw)

    # ----------------------------------------------------------------------------------------------






