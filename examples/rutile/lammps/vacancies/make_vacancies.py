
import numpy as np


num_O_vac = 500
num_Ti_vac = num_O_vac//2
num_defects = num_O_vac+num_Ti_vac

pristine_file_name = '../pristine/rutile_12x12x16.supercell'
with open(pristine_file_name,'r') as _f:
    lines = _f.readlines()

header = lines[:16]

data = lines[16:]
data = [_.strip().split() for _ in data]
data = np.array(data,dtype=float)
data = data[np.argsort(data[:,0]),:]

num_atoms = data.shape[0]


header[2] = f'{num_atoms-num_defects} atoms\n'

O_inds = np.flatnonzero(data[:,1] == 2)
Ti_inds = np.flatnonzero(data[:,1] == 1)

np.random.shuffle(O_inds)
np.random.shuffle(Ti_inds)


O_inds = O_inds[num_O_vac:]
Ti_inds = Ti_inds[num_Ti_vac:]

inds = np.r_[O_inds,Ti_inds]
inds = inds[np.argsort(inds)]

data = data[inds,:]

print(np.sum(data[:,2]))

num_atoms -= num_defects

defect_file = f'rutile_Ov_{num_O_vac:g}_Tiv_{num_Ti_vac:g}_12x12x16.supercell'
with open(defect_file,'w') as _f:
    for line in header:
        _f.write(line)
    for ii in range(num_atoms):
        _ = f'{ii+1:5g} {data[ii,1]:3g} {data[ii,2]: 5.4f} {data[ii,3]:14.9f} {data[ii,4]:14.9f} {data[ii,5]:14.9f}\n'
        _f.write(_)





