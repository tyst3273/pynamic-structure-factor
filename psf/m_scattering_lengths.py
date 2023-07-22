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
import multiprocessing as mp

from psf.m_import import import_module

class c_scattering_lengths:

    """
    holds neutron scattering lengths to look up based on atoms types and x-ray scattering form 
    factor parameters to be looked up based on atom types and later used to calculate f(Q)
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm,timers):

        """
        scattering lengths for neutron/xray scattering
        """

        self.config = config
        self.comm = comm
        self.timers = timers

        self.num_types = self.config.num_types
        self.atom_types = self.config.atom_types

        self.experiment_type = self.config.experiment_type

        # read depending on exp type
        if self.experiment_type == 'neutrons':
            self._read_neutron_scattering_lengths()
        if self.experiment_type == 'xrays':
            self._read_xray_scattering_params()

        # now go an use data for stuff
        self._get_scattering_data()

    # ----------------------------------------------------------------------------------------------

    def _get_scattering_data(self):
        
        """
        set up the scattering 'data' i.e. either scattering lengths or xray form factors
        calculate on Q-point grid if xrays since form factors depend on Q
        """

        msg = '*** experiment type ***\n'
        if self.experiment_type == 'xrays':
            msg += 'the experiment is using xrays!\n'
            print(msg)
            self._get_xray_form_factors()

        if self.experiment_type == 'neutrons':
            msg += 'the experiment is using neutrons!'
            print(msg)
            self._get_neutron_scattering_lengths()

    # ----------------------------------------------------------------------------------------------

    def _get_xray_form_factors(self):

        """
        calculate x-ray from factors on Q-points

        compute xray form factor f(|Q|) using the data from the params dictionary.
        f(|Q|) = c+sum_{i=(1,2,3,4)} ai*exp(-bi*(|Q|/4/pi)**2)

        NOTE: this is now parallelized for speed!
        I should add a low memory option that calculates only per type and map the 
        the types to atoms. the calc. will be slower but...
        """

        self.timers.start_timer('get_form_factors',units='s')

        msg = 'calculating x-ray form factors. this might take a while ...\n'
        print(msg)

        # form factor array
        _num_atoms = self.comm.traj.num_atoms
        _num_Q = self.comm.Qpoints.num_Q

        self.form_factors = np.zeros((_num_atoms,_num_Q))
        
        # get stuff for parallelization
        _num_procs = self.comm.paral.num_Qpoint_procs

        # Queue for passing data between procs
        self.queue = mp.Queue()

        # set up parallelism
        _procs = []
        for _proc in range(_num_procs):
            _procs.append(mp.Process(target=self._form_factors_on_proc,args=[_proc]))

        # start execution
        for _proc in _procs:
            _proc.start()

        # get the results from queue
        for _proc in range(_num_procs):
            _form_factors, _proc = self.queue.get()
            _Q_inds = self.comm.paral.Qpts_on_proc[_proc]
            self.form_factors[:,_Q_inds] = _form_factors[...]

        # close queue, kill it
        self.queue.close()
        self.queue.join_thread()

        # blocking; wait for all procs to finish before moving on
        for _proc in _procs:
            _proc.join()
            _proc.close()

        self.timers.stop_timer('get_form_factors')

    # ----------------------------------------------------------------------------------------------

    def _form_factors_on_proc(self,proc=0):

        """
        calculate x-ray from factors on Q-points

        compute xray form factor f(|Q|) using the data from the params dictionary.
        f(|Q|) = c+sum_{i=(1,2,3,4)} ai*exp(-bi*(|Q|/4/pi)**2)
        """

        _num_atoms = self.comm.traj.num_atoms
        _Qpt_inds = self.comm.paral.Qpts_on_proc[proc]
        _Q = self.comm.Qpoints.Q_len[_Qpt_inds]
        _num_Q = _Q.size
        _params = self.scattering_params
        _num_types = self.config.num_types
        _types = self.comm.traj.types

        # form factor for all Q on this proc
        _ffacs_Q = np.zeros((_num_atoms,_num_Q))

        # work array
        _ffacs = np.zeros(_num_types)
        for ii in range(_num_Q):

            if proc == 0:
                if ii % 100 == 0:
                    _ = f'Q[{ii:g}]'
                    msg = f'now on {_:8} out of {_num_Q:g}'
                    print(msg)

            _Q4pi = (_Q[ii]/4/np.pi)**2
            
            # get the form factors for each atom type
            for jj in range(_num_types):
                _p = _params[jj]
                _ffacs[jj] = _p['c']+ \
                _p['a1']*np.exp(-_p['b1']*_Q4pi)+ \
                _p['a2']*np.exp(-_p['b2']*_Q4pi)+ \
                _p['a3']*np.exp(-_p['b3']*_Q4pi)+ \
                _p['a4']*np.exp(-_p['b4']*_Q4pi)

            # assign form factors to all atoms;
            # would be faster to use fancy indexing
            for jj in range(_num_atoms):
                _type = _types[jj]
                _ffacs_Q[jj,ii] = _ffacs[_type]

        self.queue.put([_ffacs_Q, proc])

    # ----------------------------------------------------------------------------------------------

    def _get_neutron_scattering_lengths(self):

        """
        set up the scattering 'data' i.e. either scattering lengths or xray form factors
        calculate on Q-point grid if xrays since form factors depend on Q
        """

        _num_atoms = self.comm.traj.num_atoms
        _types = self.comm.traj.types
        _xlens = self.neutron_scattering_lengths

        self.scattering_lengths = np.zeros(_num_atoms)
        for ii in range(_num_atoms):
                
            _type = _types[ii] 
            self.scattering_lengths[ii] = _xlens[_type]

    # ----------------------------------------------------------------------------------------------

    def _read_xray_scattering_params(self):

        """
        lookup xray scattering params
        """

        self.scattering_params = [] 

        _params = import_module('psf.scattering_data.xray_scattering_params')

        msg = '\n*** xray scattering params ***\n'
        msg += ' (type)   (a1)   (b1)   (a2)   (b2)   (a3)   (b3)   (a4)   (b4)   (c)  \n'
        for ii in range(self.num_types):

            _ = {}
            
            _x = _params.scattering_params[self.atom_types[ii]]
            _['a1'] = _x['a1']
            _['b1'] = _x['b1']
            _['a2'] = _x['a2']
            _['b2'] = _x['b2']
            _['a3'] = _x['a3']
            _['b3'] = _x['b3']
            _['a4'] = _x['a4']
            _['b4'] = _x['b4']
            _['c'] = _x['c']

            self.scattering_params.append(_)

            msg += f'  {self.atom_types[ii]:4}  {_x["a1"]:6.2f} {_x["b1"]:6.2f}' \
                f' {_x["a2"]:6.2f} {_x["b2"]:6.2f} {_x["a3"]:6.2f} {_x["b3"]:6.2f}' \
                f' {_x["a4"]:6.2f} {_x["b4"]:6.2f} {_x["c"]:6.2f}\n'
        
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _read_neutron_scattering_lengths(self):

        """
        lookup neutron scattering lengths
    
        NOTE: do I actually need to discard the imaginary part? Since S(Q) ~ |b|^2, we will still
        get real data ... it will just be more computational expensive to have b complex ...
        if it comes up, I can fix it. 
        """

        eps = 0.0001
        
        self.neutron_scattering_lengths = np.zeros(self.num_types)

        _xlens = import_module('psf.scattering_data.neutron_scattering_lengths')

        msg = '\n*** neutron scattering lengths ***\n'
        msg += ' (type)  (b in fm)\n'
        for ii in range(self.num_types):
            
            _x = _xlens.scattering_lengths[self.atom_types[ii]]
            if np.abs(np.imag(_x)) > eps:
                msg += '\nWARNING! the neutron scattering lenght for type ' \
                    f'\'{self.atom_types[ii]}\' has a large\nimaginary part! ' \
                     'i will discard the imaginary part but the results\nmay not be sensible ...\n\n'
                print(msg)
            _x = np.real(_x)
            self.neutron_scattering_lengths[ii] = _x

            msg += f'  {self.atom_types[ii]:4}  {_x:8.4f}\n'  

        print(msg)

    # ----------------------------------------------------------------------------------------------


