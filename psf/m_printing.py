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

# --------------------------------------------------------------------------------------------------

def print_preamble():

    """
    self explanatory ...
    """

    preamble = '\n\n#######################################################################\n'
    preamble += 'Pynamic Structure Factor \n'
    preamble += 'author: Tyler C. Sterling\n'
    preamble += 'email: ty.sterling@colorado.edu\n'
    preamble += 'affil: Physics Dept., University of Colorado Boulder\n'
    preamble += '  Neurtron Scattering and Raman Spectroscopy Lab\n'
    preamble += '#######################################################################\n'

    print(preamble)

# --------------------------------------------------------------------------------------------------

def print_goodbye():

    """
    self explanatory ...
    """

    goodbye = '#######################################################################\n'
    goodbye += 'the calculation finished willy-nilly\n'
    goodbye += 'as always, check the results carefully\n'
    goodbye += 'bye!\n'
    goodbye += '#######################################################################\n'

    print(goodbye)

# --------------------------------------------------------------------------------------------------

