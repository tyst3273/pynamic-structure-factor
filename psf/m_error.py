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
        print('*** exception ***')
        print(str(exception)+'\n')

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





