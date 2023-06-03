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

import psf.m_lattice as m_lattice
import psf.m_Qpoints as m_Qpoints
import psf.m_timing as m_timing
import psf.m_scattering_lengths as m_scattering_lengths
import psf.m_trajectory as m_trajectory
import psf.m_structure_factors as m_structure_factors
import psf.m_parallel as m_parallel
import psf.m_io as m_io


class c_communicator:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,config,timers):
        
        """
        put other objects inside of this to conveniently pass them in/out of classes and methods
        """

        self.config = config
        self.timers = timers

    # ----------------------------------------------------------------------------------------------

    def setup_calculation(self,pos=None,types=None):

        """
        initialize errything
        """

        # lattice and reciprocal lattice vectors
        self.lattice = m_lattice.c_lattice(self.config,self)

        # Q-points created using one of numerous methods specified in config
        self.Qpoints = m_Qpoints.c_Qpoints(self.config,self,self.timers)
        self.Qpoints.generate_Qpoints()

        # set up parallelism over Q-points
        self.paral = m_parallel.c_multi_processing(self.config,self)

        # find trajectory file and set up for block-averaging
        self.traj = m_trajectory.c_trajectory(self.config,self,self.timers)

        # setup using external data if requested
        if types is not None:
            self.traj.set_external_types(types)
        if pos is not None:
            self.traj.set_external_pos(pos)

        # scattering lengths for the S(Q,w) calculation
        self.xlengths = m_scattering_lengths.c_scattering_lengths(self.config,self,self.timers)

        # setup object to hold structure factors, calculate stuff, etc.
        self.strufacs = m_structure_factors.c_structure_factors(self.config,self,self.timers)

    # ----------------------------------------------------------------------------------------------



