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
import sys
import importlib
import importlib.machinery as machinery

# custom modules
from psf.m_error import crash


# --------------------------------------------------------------------------------------------------

def import_module(module_name):

    """
    import a module that is located in the source 'tree'. 
    this doenst suffer from the issue below since these modules all have unique names!
    """

    try:
        module = importlib.import_module(module_name)
    except Exception as _ex:
        crash(f'the file \'{module_name}\' is broken\n',_ex)

    return module

# --------------------------------------------------------------------------------------------------

def import_module_from_path(module_path):

    """
    import a module with arbitrary path, e.g. ./input_params.py 
    hopefully this method will be robust...

    UPDATE: apparently modules are global. loading a new module with the same name as an
    old one just adds new data to the old one... it doesnt delete it and its attributes like
    I hoped! removing it from sys.modules seems to work.
    """

    if 'psf_module' in sys.modules:
        sys.modules.pop('psf_module')

    try:
        module = importlib.machinery.SourceFileLoader('psf_module',module_path)
        module = module.load_module()
    except Exception as ex:
        crash(f'the file \'{module_path}\' is broken\n',ex)

    return module

# --------------------------------------------------------------------------------------------------




