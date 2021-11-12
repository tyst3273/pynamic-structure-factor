
# system modules
import os
import traceback
import sys

# custom modules

class c_logger:

    """
    class to handle writing to files and standard out. probably not the cleanest organization,
    but i also put the wrapper to exit in here. call c_logger.crash([msg]) to print a message
    and traceback then exit
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self):

        """
        open the files
        """

        self.info_path = 'INFO'
        self.debug_path = 'DEBUG'
        self.crash_path = 'CRASH'
        self._info_file = open(self.info_path,'w')
        self._debug_file = open(self.debug_path,'w')
        self._crash_file = open(self.crash_path,'w')

    # ----------------------------------------------------------------------------------------------

    def console(self,msg):

        """
        unformatted print to screen 
        """
        
        self._print(msg)

    # ----------------------------------------------------------------------------------------------

    def info(self,msg):

        """
        print to info file
        """
                
        self._print(msg,msg_type='INFO')
        self._print(msg,handle=self._info_file)

    # ----------------------------------------------------------------------------------------------

    def debug(self,msg):

        """
        write debug info to debug file
        """

        self._print(msg,handle=self._debug_file)

    # ----------------------------------------------------------------------------------------------

    def crash(self,msg='an error occured'):
        
        """
        print and write error message, trace back, and crash
        """
        
        self._print(msg,msg_type='CRASH')
        self._print(msg,handle=self._crash_file)
        tb = traceback.print_stack(file=self._crash_file)
        self.clean_up()
        self._exit()

    # ----------------------------------------------------------------------------------------------

    def _print(self,msg,handle=None,msg_type=None):

        """
        wrapper to control printing
        """

        if msg_type == None:
            msg = ' '+msg
        else:
            msg = f' ** {msg_type} ** \n '+msg
        print(msg,file=handle,end='\n\n')

    # ----------------------------------------------------------------------------------------------

    def _exit(self):

        """
        wrapper to crash. may be needed for multiprocessing to prevent deadlock
        """
        
        raise SystemExit

    # ----------------------------------------------------------------------------------------------

    def clean_up(self):

        """
        close all files and remove empty ones
        """
        
        self._info_file.close()
        self._debug_file.close()
        self._crash_file.close()

        if os.path.getsize(self.debug_path) == 0:
            os.remove(self.debug_path)
        if os.path.getsize(self.info_path) == 0:
            os.remove(self.info_path)
        if os.path.getsize(self.crash_path) == 0:
            os.remove(self.crash_path)

    # ----------------------------------------------------------------------------------------------
    # other printing utilities
    # ----------------------------------------------------------------------------------------------

    def print_bar(self,symbol='-',length=72):

        msg = symbol*length
        self._print(msg)

    # ----------------------------------------------------------------------------------------------

    def print_herald(self):

        herald = '\n' \
        ' pynamic-structure-factor (psf) \n\n' \
        ' this code can be used to calcule bragg, diffuse, and inelastic neutron \n' \
        ' and xray scattering spectra from molecular dynamics trajectories. \n\n' \
        ' author:     Tyler C. Sterling \n' \
        ' affilition: University of Colorado Boulder \n' \
        ' email:      ty.sterling@colorado.edu \n\n' \
        ' license:    not yet ... \n' \
        ' version:    0.4 ' \

        self._print(herald)
        self.print_bar()

    # ----------------------------------------------------------------------------------------------


        




        

