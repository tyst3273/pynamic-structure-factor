
import psf.m_io as m_io

reader = m_io.c_reader('./out/psf_STRUFACS.hdf5')
reader.read_sqw()
reader.read_bragg()
reader.read_diffuse()
reader.read_rho_squared()





