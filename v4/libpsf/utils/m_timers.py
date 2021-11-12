
# system modules
from timeit import default_timer

# custom modules

class _timer:

    """
    private timer contained in c_timers.timers
    """
    
    # ----------------------------------------------------------------------------------------------

    def __inif__(self,label,units,log_level):
        
        if units == 'ms':
            self.units = units
            self.scale = 1e3
        if units == 'm':
            self.units = units
            self.scale = 1/60
        else:
            self.units = 's'
            self.scale = 1
        
        self.log_level = log_level
        self.label = label
        self.start_time = 0
        self.stop_time = 0
        self.total_time = 0
        self.lap_time = 0
        self.running = False
        self.string = ' - {self.label:<24} [{self.running}] {} {} ... '

    # ----------------------------------------------------------------------------------------------

    def start(self):

        self.start_time = default_timer()

    # ----------------------------------------------------------------------------------------------

    def stop(self):

        self.stop_time = default_timer()
        self.lap_time = (self.stop_time-self.start_time)*self.scale
        self.total_time = self.total_time+self.lap_time

    # ----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# the timers class used in the code
# --------------------------------------------------------------------------------------------------

class c_timers:

    """
    class to hold and run timers
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self):

        self.timers = {}

    # ----------------------------------------------------------------------------------------------
        
    def create_timer(self,label,units='s',log_level='info'):

        self.timers[label] = _timer(label,units=units,log_level=log_level)
    
    # ----------------------------------------------------------------------------------------------

    def start_timer(self,label):

        pass

    # ----------------------------------------------------------------------------------------------

    def stop_timer(self,label):

        pass

    # ----------------------------------------------------------------------------------------------

    def print_timing(self,label):

        pass

    # ----------------------------------------------------------------------------------------------


