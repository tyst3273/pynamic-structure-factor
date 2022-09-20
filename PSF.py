#!/home/ty/anaconda3/bin/python3

"""
'driver' script for psf code. gonna put 'tasks' in here, i.e. macros for calculating 
all kinds of stuff in here... but for now, its just the main script.
"""

import psf.m_config as m_config
import psf.m_lattice as m_lattice
import psf.m_Qpoints as m_Qpoints

preamble = '\n\n#######################################################################\n'
preamble += 'Pynamic Structure Factors \n'
preamble += 'author: Tyler C. Sterling\n'
preamble += 'email: ty.sterling@colorado.edu\n'
preamble += 'affil: Physics Dept., University of Colorado Boulder\n'
preamble += '  Neurtron Scattering and Raman Spectroscopy Lab\n'
preamble += '#######################################################################\n'
print(preamble)

config = m_config.c_config()
config.get_args_from_file()

lattice = m_lattice.c_lattice(config)


#Qpts = m_Qpoints.c_Qpoints(config,lattice)



