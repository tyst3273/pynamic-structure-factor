#!/home/ty/anaconda3/bin/python3

"""
'driver' script for psf code. gonna put 'tasks' in here, i.e. macros for calculating 
all kinds of stuff in here... but for now, its just the main script.
"""

import psf.m_communicator as m_communicator
import psf.m_config as m_config
import psf.m_lattice as m_lattice
import psf.m_Qpoints as m_Qpoints
import psf.m_timing as m_timing
import psf.m_scattering_lengths as m_scattering_lengths

preamble = '\n\n#######################################################################\n'
preamble += 'Pynamic Structure Factors \n'
preamble += 'author: Tyler C. Sterling\n'
preamble += 'email: ty.sterling@colorado.edu\n'
preamble += 'affil: Physics Dept., University of Colorado Boulder\n'
preamble += '  Neurtron Scattering and Raman Spectroscopy Lab\n'
preamble += '#######################################################################\n'
print(preamble)



comm = m_communicator.c_communicator()

timers = m_timing.c_timers()
timers.start_timer('PSF',units='m')

config = m_config.c_config()
config.get_args_from_file()


comm.lattice = m_lattice.c_lattice(config,comm)


comm.Qpoints = m_Qpoints.c_Qpoints(config,comm,timers)
comm.Qpoints.generate_Qpoints()


comm.xlengths = m_scattering_lengths.c_scattering_lengths(config,comm)


timers.stop_timer('PSF')
timers.print_timing()





