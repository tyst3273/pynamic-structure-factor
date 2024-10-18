import h5py

with h5py.File('pos.hdf5','r') as db:
    box_vecs = db['box_vectors'][...]

lat_vecs = box_vecs.mean(axis=0)
lat_vecs[0,:] /= 10
lat_vecs[1,:] /= 10
lat_vecs[2,:] /= 12
print(lat_vecs)


