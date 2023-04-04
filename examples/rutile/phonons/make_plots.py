
from psf.m_io import c_reader

reader = c_reader('./out/pristine_neutrons_STRUFACS.hdf5')
reader.read_sqw()
Q_rlu = reader.Q_rlu
pn_sqw = reader.sqw

reader = c_reader('./out/pristine_neutrons_STRUFACS.hdf5')
reader.read_sqw()
Q_rlu = reader.Q_rlu
sn_sqw = reader.sqw




