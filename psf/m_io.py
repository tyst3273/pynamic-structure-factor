
import os
import numpy as np
import h5py 


class c_io:
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,comm,timers):

        """
        read/write file containing structure factors
        """

        # copy refs to stuff
        self.comm = comm
        self.config = config

        # whether or no 'header' has been written
        self.header_flag = False

        # get the file name
        self.output_prefix = self.config.output_prefix
        self.output_directory = self.config.output_directory

        print(self.output_prefix)

    # ----------------------------------------------------------------------------------------------

    def _write_header(self):

        """
        write S(Q,w) info to file
        """

        if not self.header_flag:
            self._write_header()
            self.header_flag = True

        pass

    # ----------------------------------------------------------------------------------------------
    
    def write_sqw(self):

        """
        write S(Q,w) info to file
        """

        if not self.header_flag:
            self._write_header()
            self.header_flag = True

        pass

    # ----------------------------------------------------------------------------------------------

    def read_sqw(self):

        """
        read S(Q,w) info from file
        """

        pass

    # ----------------------------------------------------------------------------------------------

    def write_bragg(self):

        """
        write bragg info to file
        """

        if not self.header_flag:
            self._write_header()
            self.header_flag = True

        pass

    # ----------------------------------------------------------------------------------------------

    def read_bragg(self):

        """
        read bragg info from file
        """

        pass

    # ----------------------------------------------------------------------------------------------

    def write_diffuse(self):

        """
        write S(Q,w) info to file
        """

        if not self.header_flag:
            self._write_header()
            self.header_flag = True

        pass

    # ----------------------------------------------------------------------------------------------

    def read_diffuse(self):

        """
        read S(Q,w) info from file
        """

        pass

    # ----------------------------------------------------------------------------------------------

