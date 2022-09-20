
# system modules
import os

# custom modules


# --------------------------------------------------------------------------------------------------

def kill():

    """
    kill the program
    """

    raise SystemExit

# --------------------------------------------------------------------------------------------------

def crash(msg=None,exception=None):

    """
    kill the program
    """
        
    print('\n*** error ***')

    if not msg is None:
        print(msg)

    if not exception is None:
        print(exception)

    kill()

# --------------------------------------------------------------------------------------------------

def check_file(file_path,warn=False):

    """
    check if a file is missing
    """

    if not os.path.exists(file_path):
        msg = f'the file:\n  \'{file_path}\'\nis missing'
        if warn:
            print('\n*** warning! ***\n'+msg+'\n\ncontinuing but might crash!\n')
        else:
            crash(msg+'\n')

# --------------------------------------------------------------------------------------------------





