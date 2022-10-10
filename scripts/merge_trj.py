

from psf.utils.m_mergers import c_user_data

num_atoms = 49152

for T in [200,230,250,300,350]:

    trj = c_user_data(number_of_atoms=num_atoms,output_file=f'mapi_{T:g}K.hdf5',
        input_files=f'dump_MAPbI3_16x16x16_T{T:g}.xyz',file_format='lammpstrj',block_size=10)
    trj.merge_files()







