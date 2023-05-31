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

# system modules
from timeit import default_timer

# custom modules
from psf.m_error import crash


class _timer:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,label,units):

        """
        private class for the timers class to use. holds timer objects with all 
        their stuff set in them
        """
        if units == 'ms':
            self.scale = 1e3
            self.units = 'ms'
        elif units == 'm':
            self.scale = 1/60
            self.units = 'm'
        else:
            self.scale = 1
            self.units = 's'
        self.label = label
        self.start_time = default_timer()
        self.lap_time = 0
        self.total_time = 0
        self.running = False
        self.calls = 0
        self.set_str()

    # ----------------------------------------------------------------------------------------------

    def set_str(self):
        
        """
        create string for timer to print
        """

        self.str = f'  {self.label:<24} {self.calls:8g} {self.lap_time: 10.3f} ' \
               f'{self.total_time: 10.3f}   {self.units:<2}'

    # ----------------------------------------------------------------------------------------------

    def start(self):
        
        """
        start the timer
        """

        if self.running:
            self.stop()
        self.start_time = default_timer()
        self.running = True

    # ----------------------------------------------------------------------------------------------

    def stop(self):

        """
        stop the timer
        """

        self.lap_time = (default_timer()-self.start_time)*self.scale
        self.total_time = self.total_time+self.lap_time
        self.running = False
        self.calls += 1
        self.set_str()
        
    # ----------------------------------------------------------------------------------------------

    def print_timing(self):

        """
        print timing info if timer used locally
        """

        if self.running: 
            msg = f'WARNING! the timer: \'{self.label}\' ' \
                          'is still running!'
            print(msg)

        msg = '\nlocal timer:\n'
        msg += '  (label)                    (calls)     (lap)     (total) \n'
        msg += self.str+'\n'
        print(msg)

    # ----------------------------------------------------------------------------------------------


class c_timers:

    # ----------------------------------------------------------------------------------------------

    def __init__(self):

        """
        create timers for whatever
        """

        self.timers = {}

    # ----------------------------------------------------------------------------------------------

    def start_timer(self,label,units='s'):

        """
        start the requested timer. create it if it doesnt exist
        """
    
        # create if it doesnt exists
        if label not in self.timers.keys():
            self.timers[label] = _timer(label,units)

        # start it
        self.timers[label].start()

    # ----------------------------------------------------------------------------------------------

    def stop_timer(self,label):

        """
        stop the requested timer
        """

        self.timers[label].stop()

    # ----------------------------------------------------------------------------------------------

    def stop_all_running_timers(self):

        """
        stop all running timers
        """

        for label in self.timers.keys():
            if self.timers[label].running:
                self.timers[label].stop()

    # ----------------------------------------------------------------------------------------------

    def print_timing(self):

        """
        print all the timing info to the info file
        """

        timing = '\n*** timing *** \n' \
              '  (label)                    (calls)     (lap)     (total) \n'
        for label in self.timers.keys():
            if self.timers[label].running:
                timing += f'WARNING! the timer: \'{self.timers[label].label}\' ' \
                          'was still running!\n'
                self.timers[label].stop()                
            timing += self.timers[label].str+'\n'
        print(timing)

    # ----------------------------------------------------------------------------------------------




