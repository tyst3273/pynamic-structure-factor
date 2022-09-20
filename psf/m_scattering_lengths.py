
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

        # get data depending on exp type
        if self.config.experiment_type == 'neutrons':
            self._get_neutron_scattering_lengths()
        if self.config.experiment_type == 'xrays':
            self._get_xray_scattering_params()
  
    # ----------------------------------------------------------------------------------------------

    def _get_xray_scattering_params(self):

        """
        lookup xray scattering params
        """

        self.scattering_params = np.zeros((self.num_types,9)) # 9 params per type

        _params = import_module('psf.scattering_data.xray_scattering_params')

        msg = '\n*** xray scattering params ***\n'
        msg += ' (type)   (a1)   (b1)   (a2)   (b2)   (a3)   (b3)   (a4)   (b4)   (c)  \n'
        for ii in range(self.num_types):
            
            _x = _params.scattering_params[self.unique_types[ii]]
            self.scattering_params[ii,0] = _x['a1']
            self.scattering_params[ii,1] = _x['b1']
            self.scattering_params[ii,2] = _x['a2']
            self.scattering_params[ii,3] = _x['b2']
            self.scattering_params[ii,4] = _x['a3']
            self.scattering_params[ii,5] = _x['b3']
            self.scattering_params[ii,6] = _x['a4']
            self.scattering_params[ii,7] = _x['b4']
            self.scattering_params[ii,8] = _x['c']

            msg += f'  {self.unique_types[ii]:4}  {_x["a1"]:6.2f} {_x["b1"]:6.2f}' \
                f' {_x["a2"]:6.2f} {_x["b2"]:6.2f} {_x["a3"]:6.2f} {_x["b3"]:6.2f}' \
                f' {_x["a4"]:6.2f} {_x["b4"]:6.2f} {_x["c"]:6.2f}\n'
        
        print(msg)


    # ----------------------------------------------------------------------------------------------

    def _get_neutron_scattering_lengths(self):

        """
        lookup neutron scattering lengths
        """

        eps = 0.0001
        
        self.scattering_lengths = np.zeros(self.num_types)

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
            self.scattering_lengths[ii] = _x

            msg += f'  {self.unique_types[ii]:4}  {_x:8.4f}\n'  

        print(msg)

    # ----------------------------------------------------------------------------------------------


