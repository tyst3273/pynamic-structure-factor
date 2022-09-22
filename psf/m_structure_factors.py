
import numpy as np


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
        self.calc_timeavg = config.calc_timeavg
        self.calc_bragg = config.calc_bragg

        # set stuff up
        if self.calc_sqw:
            self._setup_energies()
            self.sqw = np.zeros((self.comm.Qpoints.num_Q,self.num_energy),dtype=float)

            msg = '\n*** dynamic structure factor ***\ncalculating S(Q,E)...\n'
            msg += f'num_energy: {self.num_energy}\n'
            msg += f'energy_step:\n {self.energy_step: 8.6f} (meV)\n'
            msg += f' {self.energy_step/_thz2meV: 8.6f} (tHz)\n'
            msg += f'energy_max:\n {self.energy_max: 8.4f} (meV)\n'
            msg += f' {self.energy_max/_thz2meV: 8.4f} (tHz)'
            print(msg)

        if self.calc_timeavg:
            self.sq_avg = np.zeros(self.comm.Qpoints.num_Q,dtype=float)
            msg = '\n*** time-avg structure factors ***\ncalculating diffuse S(Q)\n'
            msg += '  ~ <|F(Q,t)|**2> '
            print(msg)

        if self.calc_bragg:
            self.sq_bragg = np.zeros(self.comm.Qpoints.num_Q,dtype=float)
            msg = '\n*** bragg structure factors ***\ncalculating bragg intensity S(Q)\n'
            msg += '  ~ |<F(Q),t)>|**2 '
            print(msg)

    # ----------------------------------------------------------------------------------------------

    def _setup_energies(self):

        """
        get energies based on size of blocks of traj
        """

        self.num_energy = self.comm.traj.num_steps

        # timestep is in fs, convert to freq. in tHz, then to meV
        self.energy = np.fft.fftfreq(self.num_energy,self.comm.traj.md_time_step*1e-3)*_thz2meV
        self.energy_step = self.energy[1]-self.energy[0]
        self.energy_max = self.energy.max()

    # ----------------------------------------------------------------------------------------------

    

