#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                                           !
#   ! Copyright 2025 by Tyler C. Sterling                                       !
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


