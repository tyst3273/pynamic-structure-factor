
# system modules
import argparse
import configparser
import os
import numpy as np

# custom modules

class c_config:

    """
    get command line and configuration file args
    """
    
    # ------------------------------------------------------------------------------------------

    def __init__(self,logger):
        
        """
        get command line args and set default config options

        to add a new arg, add its handle to self.kw and add it to the arg_parser. it has to be 
        added to kw because we preprocess it and append '--' to all kw. otherwise argparser is
        unhappy.
        """

        self.logger = logger

        # ------------------------------------------------------------------------------------------
        # set and get command line args
        # ------------------------------------------------------------------------------------------

        cmd_parser = argparse.ArgumentParser(description='documentation for psf command line args')

        cmd_parser.add_argument('-i','--input_file',default=['psf.ini'],
                            help='input filename(s) for psf.py',nargs='+')
        cmd_parser.add_argument('-n','-np','--num_processes',default=1,type=int,
                            help='number of \'multiprocessing\' processes to use')

        cmd_args = cmd_parser.parse_args()

        self.input_files = cmd_args.input_file
        self.num_files = len(self.input_files)
        self.num_process = cmd_args.num_processes

        # check if input files exist
        crash = False
        msg = 'the following files are missing'
        for in_file in self.input_files:
            if not os.path.exists(in_file):
                crash = True
                msg = msg+f'\n  - {in_file:<}' 
        if crash:
            self.logger.crash(msg)

        # print the files to do
        msg = f'there are {self.num_files} files to do'
        for in_file in self.input_files:
            msg = msg+f'\n  - {in_file:<}'       
            self.logger.console(msg)

        # ------------------------------------------------------------------------------------------
        # set config file opts. gotten later
        # ------------------------------------------------------------------------------------------

        # put all args in this list. order doesnt matter but theyre needed
        self.kw = ['lattice_vectors','lattice_scale','box_vectors','box_scale']

        # now set up config file defaults. set a default for everything
        self.arg_parser = argparse.ArgumentParser('documentation for the input variables')

        help_msg = 'lattice vectors of the crystal in angstrom. \n' \
                    'should be given as row vectors'
        self.arg_parser.add_argument('--lattice_vectors',nargs=9,type=float,help=help_msg,
                            default=[1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0])

        help_msg = 'homogeneous scaling for the lattice vectors' 
        self.arg_parser.add_argument('--lattice_scale',default=5.431,nargs=1,type=float,
                            help=help_msg)

        help_msg = 'simulation box vectors. transformation from primitive cell vectors \n' \
                   'is calculated from this.' 
        self.arg_parser.add_argument('--box_vectors',nargs=9,type=float,help=help_msg,
                            default=[1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0])

        help_msg = 'homogeneous scaling for the box vectors' 
        self.arg_parser.add_argument('--box_scale',default=5.431,nargs=1,type=float,
                            help=help_msg)

    # ----------------------------------------------------------------------------------------------

    def _set_args(self,trimmed_line):

        """
        use argparser to get args from string
        """

        args = self.arg_parser.parse_args(trimmed_line)
        self.lattice_vectors = np.array(args.lattice_vectors).reshape((3,3))

    # ----------------------------------------------------------------------------------------------

    def get_args(self,input_file):

        """
        read the whole text file in as one string. this thing is agnostic of the variable name
        """
        
        # check if input file can be read
        try:
            with open(input_file,'r') as f_in:
                read_lines = f_in.readlines()
        except:
            msg = f'the file \'{input_file}\' is broken'
            self.logger.crash(msg)

        # remove blank and comment lines
        lines = []
        for line in read_lines:
            if len(line.split()) == 0 or line.strip().startswith('#'):
                continue
            else:
                line = line.split('#')[0].strip().split()
                lines.extend(line)
        
        # remove '=' and add a '--' to each kw so argparser won't bitch about it
        trimmed_line = []
        for arg in lines:
            arg = arg.strip('\'" \\n ')
            if arg in ['=',':']:
                continue
            elif arg.lower() in self.kw:
                arg = '--'+arg
            trimmed_line.append(arg)

        # now go and write the args from the arg parser to attributes of this class
        self._set_args(trimmed_line)

    # ----------------------------------------------------------------------------------------------


