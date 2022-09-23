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
from psf.m_import import import_module

class c_scattering_lengths:

    """
    holds neutron scattering lengths to look up based on atoms types and x-ray scattering form 
    factor parameters to be looked up based on atom types and later used to calculate f(Q)
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm):

        """
        scattering lengths for neutron/xray scattering
        """

        self.config = config
        self.comm = comm

        self.num_types = self.config.num_types
        self.unique_types = self.config.unique_types

        self.experiment_type = self.config.experiment_type

        # read depending on exp type
        if self.experiment_type == 'neutrons':
            self._read_neutron_scattering_lengths()
        if self.experiment_type == 'xrays':
            self._read_xray_scattering_params()

        # now go an use data for stuff
        self.get_scattering_data()

    # ----------------------------------------------------------------------------------------------

    def get_scattering_data(self):
        
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
        """

        _num_atoms = self.comm.traj.num_atoms
        _num_Q = self.comm.Qpoints.num_Q
        _Q = self.comm.Qpoints.Q_len
        _params = self.scattering_params
        _num_types = self.config.num_types
        _types = self.comm.traj.types

        self.form_factors = np.zeros((_num_atoms,_num_Q))

        _form_facs = np.zeros(_num_types)
        for ii in range(_num_Q):

            _Q4pi = (_Q[ii]/4/np.pi)**2
            
            # get the form factors for each atom type
            for jj in range(_num_types):
                _form_facs[jj] = _params[jj]['c']+ \
                _params[jj]['a1']*np.exp(-_params[jj]['b1']*_Q4pi)+ \
                _params[jj]['a2']*np.exp(-_params[jj]['b2']*_Q4pi)+ \
                _params[jj]['a3']*np.exp(-_params[jj]['b3']*_Q4pi)+ \
                _params[jj]['a4']*np.exp(-_params[jj]['b4']*_Q4pi)

            # assign form factors to all atoms
            for jj in range(_num_atoms):
                _type = _types[jj]
                self.form_factors[jj,ii] = _form_facs[_type]

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

        self.scattering_params = [] #np.zeros((self.num_types,9)) # 9 params per type

        _params = import_module('psf.scattering_data.xray_scattering_params')

        msg = '\n*** xray scattering params ***\n'
        msg += ' (type)   (a1)   (b1)   (a2)   (b2)   (a3)   (b3)   (a4)   (b4)   (c)  \n'
        for ii in range(self.num_types):

            _ = {}
            
            _x = _params.scattering_params[self.unique_types[ii]]
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

            msg += f'  {self.unique_types[ii]:4}  {_x["a1"]:6.2f} {_x["b1"]:6.2f}' \
                f' {_x["a2"]:6.2f} {_x["b2"]:6.2f} {_x["a3"]:6.2f} {_x["b3"]:6.2f}' \
                f' {_x["a4"]:6.2f} {_x["b4"]:6.2f} {_x["c"]:6.2f}\n'
        
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def _read_neutron_scattering_lengths(self):

        """
        lookup neutron scattering lengths
        """

        eps = 0.0001
        
        self.neutron_scattering_lengths = np.zeros(self.num_types)

        _xlens = import_module('psf.scattering_data.neutron_scattering_lengths')

        msg = '\n*** neutron scattering lengths ***\n'
        msg += ' (type)  (b in fm)\n'
        for ii in range(self.num_types):
            
            _x = _xlens.scattering_lengths[self.unique_types[ii]]
            if np.abs(np.imag(_x)) > eps:
                msg += '\nWARNING! the neutron scattering lenght for type ' \
                    f'\'{self.unique_types[ii]}\' has a large\nimaginary part! ' \
                     'i will discard the imaginary part but the results\nmay not be sensible ...\n\n'
                print(msg)
            _x = np.real(_x)
            self.neutron_scattering_lengths[ii] = _x

            msg += f'  {self.unique_types[ii]:4}  {_x:8.4f}\n'  

        print(msg)

    # ----------------------------------------------------------------------------------------------


