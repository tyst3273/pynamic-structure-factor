
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

        # count the warnings
        self.num_warnings = 0

        # open the files
        self.info_path = 'INFO'
        self.debug_path = 'DEBUG'
        self.crash_path = 'CRASH'
        self.warn_path = 'WARN'
        self._info_file = open(self.info_path,'w')
        self._debug_file = open(self.debug_path,'w')
        self._crash_file = open(self.crash_path,'w')
        self._warn_file = open(self.warn_path,'w')

    # ----------------------------------------------------------------------------------------------

    def console(self,msg,msg_type=None):

        """
        print to screen with an optional message type
        """
        
        self._print(msg,msg_type=msg_type)

    # ----------------------------------------------------------------------------------------------

    def info(self,msg,msg_type=None):

        """
        print to info file. msg_type is optional and if given, only prints the type to console.
        msg in info_file is unformatted
        """
                
        self.console(msg,msg_type)
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
        
        # write to console
        self.console(msg,'CRASH')
        self.console(f'traceback info is in \'{self.crash_path}\'')

        # write to crash file
        self._print(msg,handle=self._crash_file)
        traceback.print_stack(file=self._crash_file)

        # clean up and crash
        self.clean_up()
        self._exit()

    # ----------------------------------------------------------------------------------------------

    def warn(self,msg='something is wrong',trace=False):

        """
        print and write warning msg
        """

        # count the warnings
        self.num_warnings = self.num_warnings+1 

        # inform the user
        #self.console('something unexpected has happened',msg_type='WARNING')

        # write to warning file
        self._print(msg,handle=self._warn_file)

        # optionally write traceback info too
        if trace:
            traceback.print_stack(file=self._warn_file)
            self._print(msg='',handle=self._warn_file)

    # ----------------------------------------------------------------------------------------------

    def _print(self,msg,handle=None,msg_type=None):

        """
        wrapper to control printing
        """

        if msg_type == None:
            msg = ' '+msg
        else:
            msg = f' ** {msg_type} ** \n '+msg
        print(msg,file=handle,end='\n\n',flush=True)

    # ----------------------------------------------------------------------------------------------

    def _exit(self):

        """
        wrapper to crash. may be needed for multiprocessing to prevent deadlock or to 
        save files etc when crashing. note, there is another module called 'atexit' which
        runs callbacks if the code exits
        """
        
        raise SystemExit

    # ----------------------------------------------------------------------------------------------

    def clean_up(self):

        """
        close all files and remove empty ones
        """
        
        # close the files
        self._info_file.close()
        self._debug_file.close()
        self._crash_file.close()
        self._warn_file.close()

        # delete empty files
        if os.path.getsize(self.debug_path) == 0:
            os.remove(self.debug_path)
        if os.path.getsize(self.info_path) == 0:
            os.remove(self.info_path)
        if os.path.getsize(self.crash_path) == 0:
            os.remove(self.crash_path)
        if os.path.getsize(self.warn_path) == 0:
            os.remove(self.warn_path)

    # ----------------------------------------------------------------------------------------------
    # other printing utilities
    # ----------------------------------------------------------------------------------------------

    def print_bar(self,handle=None,symbol='-',length=72):

        """
        write a 'bar' like '------------------'
        you can control the symbol and length
        """

        msg = symbol*length
        self._print(msg,handle)

    # ----------------------------------------------------------------------------------------------

    def print_herald(self):

        """
        print title and author info to the screen
        """

        herald = '\n' \
        ' pynamic-structure-factor (psf) \n' \
        ' a python code to calculate neutron and xray scattering structure factors \n' \
        ' from molecular dynamics trajectories \n\n' \
        ' please cite: \n\n' \
        ' author:     Tyler C. Sterling \n' \
        ' affilition: University of Colorado Boulder \n' \
        ' email:      ty.sterling@colorado.edu \n' \
        ' url:        https://github.com/tyst3273/pynamic-structure-factor \n' \
        ' version:    0.4 \n' \
        ' license:    not yet ... ' 

        self._print(herald)
        self._print(herald,handle=self._info_file)
        self.print_bar()

    # ----------------------------------------------------------------------------------------------

    def print_goodbye(self):

        """
        tell the user the code is exiting normally but let them know about warnings
        """

        goodbye = '' \
        'the calculation has finished willy-nilly! \n\n' \
        f' {self.num_warnings} warnings were written to \'{self.warn_path}\' \n' \
        ' as always, check the results carefully before doing any physics.' 

        self.print_bar()
        self._print(goodbye)

    # ----------------------------------------------------------------------------------------------




        

