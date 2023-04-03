
import numpy as np


defect_concentration = 0.10 # 10%
defect_mass = 72.64 # mass of Ge in amu


pristine_file_name = '../pristine/si_12x12x12.supercell'
with open(pristine_file_name,'r') as _f:
    lines = _f.readlines()

header = lines[:15]
header[3] = '2 atom types\n'
header.insert(12,f'2 {defect_mass:5.3f}\n')

data = lines[15:]
data = [_.strip().split() for _ in data]
data = np.array(data,dtype=float)
data = data[np.argsort(data[:,0]),:]

num_atoms = data.shape[0]
defect_inds = np.arange(num_atoms)
np.random.shuffle(defect_inds)
defect_inds = defect_inds[:int(defect_concentration*num_atoms)]

data[defect_inds,1] = 2


defect_file = f'si_ge_{defect_concentration:.2f}_12x12x12.supercell'
with open(defect_file,'w') as _f:
    for line in header:
        _f.write(line)
    for ii in range(num_atoms):
        _ = f'{data[ii,0]:5g} {data[ii,1]:3g} {data[ii,2]:12.9f} {data[ii,3]:12.9f} {data[ii,4]:12.9f}\n'
        _f.write(_)





