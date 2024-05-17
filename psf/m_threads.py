
# system modules
import os

# custom modules
from psf.m_command_line import get_num_threads

# --------------------------------------------------------------------------------------------------

def set_num_threads():

    """
    set the number of threads used on backend. set environement variable OMP_NUM_THREADS and, 
    if built against mkl, set mkl_num_threads. 
    """

    # user requested value
    num_threads = get_num_threads()
    os.environ['OMP_NUM_THREADS'] = str(num_threads)
    os.environ['OPENBLAS_NUM_THREADS'] = str(num_threads)
     
    try:
        
        # try to set num threads in mkl
        import mkl
        mkl.set_num_threads(num_threads)

    except Exception as ex:

        pass

        #msg = '\n*** WARNING ***\n'
        #msg += 'trying to set number of threads used by intel mkl failed.\n' 
        #print(msg)
        #print('*** exception ***\n',str(ex))

# --------------------------------------------------------------------------------------------------


