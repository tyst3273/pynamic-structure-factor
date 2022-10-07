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
from scipy.fft import fft, ifft
import multiprocessing as mp
from psf.m_timing import _timer
from psf.m_error import crash


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
        self.calc_rho_squared = config.calc_rho_squared

        # set stuff up
        if self.calc_sqw:
            self._setup_energies()
            self.sqw = np.zeros((self.comm.Qpoints.num_Q,self.num_energy),dtype=float)

            msg = '\n*** dynamic structure factor ***\n'
            msg += f'num_energy: {self.num_energy}\n'
            msg += f'energy_step:\n {self.energy_step: 8.6f} (meV)\n'
            msg += f' {self.energy_step/_thz2meV: 8.6f} (tHz)\n'
            msg += f'energy_max:\n {self.energy_max: 8.4f} (meV)\n'
            msg += f' {self.energy_max/_thz2meV: 8.4f} (tHz)'
            print(msg)

        if self.calc_rho_squared:
            self.rho_sq = np.zeros((self.comm.Qpoints.num_Q, \
                        self.comm.traj.num_block_steps),dtype=float)

            msg = '\n*** rho squared ***\ncalculating \'rho squared\'\n'
            msg += '  |exp(iQ.r(t))|**2'
            print(msg)

        if self.calc_diffuse:   
            self.sq_diffuse = np.zeros(self.comm.Qpoints.num_Q,dtype=float)

            msg = '\n*** diffuse intensity ***\ncalculating diffuse intensity \n'
            msg += '  <|exp(iQ.r(t))|**2> '
            print(msg)

        if self.calc_bragg:
            self.sq_bragg = np.zeros(self.comm.Qpoints.num_Q,dtype=float)

            msg = '\n*** bragg intensity ***\ncalculating bragg intensity\n'
            msg += '  |<exp(iQ.r(t))>|**2 '
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
            self.sq_bragg /= _num_blocks*_num_atoms
        if self.calc_diffuse:
            self.sq_diffuse /= _num_blocks*_num_atoms 
        if self.calc_sqw:
            self.sqw /= _num_blocks*_num_atoms 
        if self.calc_rho_squared:
            self.rho_sq /= _num_blocks*_num_atoms 

        self.timers.stop_timer('calc_strufacs')

    # ----------------------------------------------------------------------------------------------

    def _get_res_list(self,proc):

        """
        parse the 'results' list returned by _proc_loop_on_Q(); put the data in appropriate
        arrays
        """

        res = [proc]
        if self.calc_bragg:
            res.append(self._sq_bragg)
        if self.calc_rho_squared:
            res.append(self._rho_sq)
        if self.calc_diffuse:
            res.append(self._sq_diffuse)
        if self.calc_sqw:
            res.append(self._sqw)

        return res

    # ----------------------------------------------------------------------------------------------

    def _get_arrays_from_res(self,res):

        """
        parse the 'results' list returned by _proc_loop_on_Q(); put the data in appropriate
        arrays
        """

        # get the proc
        proc = res[0]
        res.pop(0)

        # get the Q-point indices this proc did
        Q_inds = self.comm.paral.Qpts_on_proc[proc]

        # get the data and put on correct Qpts
        if self.calc_bragg:
            self.sq_bragg[Q_inds] += res[0]
            res.pop(0)
        if self.calc_rho_squared:
            self.rho_sq[Q_inds,:] += res[0]
            res.pop(0)
        if self.calc_diffuse:
            self.sq_diffuse[Q_inds] += res[0]
            res.pop(0)
        if self.calc_sqw:
            self.sqw[Q_inds,:] += res[0]
            res.pop(0)

    # ----------------------------------------------------------------------------------------------

    def _calculate_on_block(self,_block):

        """
        reads data for a block, distributes data over procs, runs them, collects data, 
        writes it to appropriate arrays
        """

        self.timers.start_timer('calc_on_block',units='m')

        # get stuff for parallelization
        _num_procs = self.comm.paral.num_Qpoint_procs

        # get the trajectory from the file
        self.comm.traj.read_trajectory_block(_block)

        # Queue for passing data between procs
        self.queue = mp.Queue()

        _procs = []
        for _proc in range(_num_procs):
            _procs.append(mp.Process(target=self._proc_loop_on_Q,args=[_proc]))

        # start execution
        for _proc in _procs:
            _proc.start()

        # get the results from queue
        for _proc in range(_num_procs):
            _res = self.queue.get()
            self._get_arrays_from_res(_res)

        # close queue, kill it
        self.queue.close()
        self.queue.join_thread()

        # blocking; wait for all procs to finish before moving on
        for _proc in _procs:
            _proc.join()
            _proc.close()

        self.timers.stop_timer('calc_on_block')

    # ----------------------------------------------------------------------------------------------

    def _get_empty_strufac_arrays(self,_nQ,proc=0):

        """
        get empty arrays that are used repeatedly to keep from having to allocate/reallocate
        """

        _num_steps = self.comm.traj.num_block_steps
        
        if self.calc_sqw:
            _num_energy = self.num_energy
            self._sqw = np.zeros((_nQ,_num_energy),dtype=float)

        if self.calc_diffuse or self.calc_rho_squared:
            self._rho_sq = np.zeros((_nQ,_num_steps),dtype=float)

        if self.calc_diffuse:
            self._sq_diffuse = np.zeros(_nQ,dtype=float)

        if self.calc_bragg:
            self._sq_bragg = np.zeros(_nQ,dtype=float)

    # ----------------------------------------------------------------------------------------------

    def _calc_strufacs(self,exp_iQr,ind):

        """
        calculate whatever was requested ... sorry for the spaghetti logic
        """

        if self.calc_bragg:
            self._sq_bragg[ind] = np.abs((exp_iQr).mean())**2

        if self.calc_diffuse or self.calc_rho_squared:
            self._rho_sq[ind,:] = np.abs(exp_iQr)**2

        if self.calc_diffuse:
            self._sq_diffuse[ind] = self._rho_sq[ind,:].mean()

        if self.calc_sqw:
            self._sqw[ind,:] = np.abs(fft(exp_iQr))**2

    # ----------------------------------------------------------------------------------------------

    def _proc_loop_on_Q(self,proc=0):

        """
        loops over the Q-points assinged to the proc, calcs S(Q,E) etc, then puts data in queue
        """

        # convenience refs for below
        _exp_type = self.comm.xlengths.experiment_type
        _num_steps = self.comm.traj.num_block_steps
        _num_atoms = self.comm.traj.num_atoms

        # get the Q-points that this proc is supposed do
        _Qpt_inds = self.comm.paral.Qpts_on_proc[proc]
        _Qpts = self.comm.Qpoints.Q_cart[_Qpt_inds,:]
        _nQ = _Qpts.shape[0]

        # get empty arrays
        self._get_empty_strufac_arrays(_nQ)

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

            # depends on Q for xrays, but calculate earlier
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

            # get different arrays depending on what is requested
            self._calc_strufacs(_exp_iQr,ii)

        # return list of results
        res = self._get_res_list(proc)

        # put results in queue to be passed to main proc
        self.queue.put(res)
        
    # ---------------------------------------------------------------------------------------------- 

    def put_on_mesh(self):

        """
        convenience methods that put S(Q,w) on unordered list of Q onto cartesian grid of Q-points
        for symmetry reduced Q-points, it also unfolds into full reciprocal space

        this calls methods from c_Qpoints
        """

        _use_mesh = self.comm.Qpoints.use_mesh
        if not _use_mesh:
            msg = 'a Qpoint-mesh was not requested, i.e. cannot be unfolded!\n' \
                  'use \'Qpoints_option\' = \'mesh\' or \'write_mesh\' \n'
            crash(msg)

        if self.calc_sqw:
            self.sqw = self.comm.Qpoints.unfold_onto_Q_mesh(self.sqw)
        if self.calc_diffuse:
            self.sq_diffuse = self.comm.Qpoints.unfold_onto_Q_mesh(self.sq_diffuse)
        if self.calc_rho_squared:
            self.rho_sq = self.comm.Qpoints.unfold_onto_Q_mesh(self.rho_sq)
        if self.calc_bragg:
            self.sq_bragg = self.comm.Qpoints.unfold_onto_Q_mesh(self.sq_bragg)

    # ----------------------------------------------------------------------------------------------






