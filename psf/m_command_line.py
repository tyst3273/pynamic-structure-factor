
"""
NOTE: set_num_threads must be called before anything is imported and it calls get_num_threads
from here. so if any other imports have to be added here, get_num_threads must be moved.
"""

# system modules
import argparse

# --------------------------------------------------------------------------------------------------

def _parse_command_line_args():

    """
    get command line args
    """

    # get cmd line args
    description = 'command line args for \'PSF.py\''
    cmd_parser = argparse.ArgumentParser(description=description)

    # input file
    help_msg = 'input file for \'PSF.py\'. it is a python file that is imported as a module\n'
    help_msg += 'to batch run input files, give multiple in the order you want to run them'
    cmd_parser.add_argument('-i','--input-files',default=[None],help=help_msg,nargs='+')

    # number of OMP threads
    help_msg = 'number of threads to use for threaded libraries (e.g. scipy)\n'
    help_msg += 'this option sets num_mkl_threads or environment variable OMP_NUM_THREADS\n'
    help_msg += 'before importing anything so that the right number of threads are used.\n'
    help_msg += 'works with multiprocessing! default is 1 (no OMP).'
    cmd_parser.add_argument('-n','--num-threads',default=1,help=help_msg,nargs='?',type=int)

    # get cmd line args
    cmd_args = cmd_parser.parse_args()

    input_files = cmd_args.input_files
    num_threads = cmd_args.num_threads

    return input_files, num_threads

# --------------------------------------------------------------------------------------------------

def get_input_files():

    """
    get input file(s) from command line args. default is None (dont use).
    """

    input_files, _ = _parse_command_line_args()
    return input_files

# --------------------------------------------------------------------------------------------------

def get_num_threads():

    """
    get number of OMP threads from command line. default is 1
    """

    _, num_threads = _parse_command_line_args()
    return num_threads

# --------------------------------------------------------------------------------------------------


