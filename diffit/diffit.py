
import numpy as np

from PSF import c_PSF
import mods.m_crystals as m_crystals

def get_file_str(ind):
    return f'_{ind}.hdf5'

# list 'c_rutile' objects that model instances of rutile
supercells = []

# target supercell size and steps around it to remove finite size effects
num_reps_target = [15,15,15]
d_reps = 5


# set up supercells from 'c_rutile' objects
for _r in range(-d_reps,d_reps+1):
    
    _reps = np.array(num_reps_target)+_r

    _rutile = m_crystals.c_rutile()
    _rutile.make_supercell(reps=_reps)

    supercells.append(_rutile)



# go and calc scattering intensity from each supercell
for ind, sc in enumerate(supercells):
    
    # object to calculate scattering intensity each time
    PSF = c_PSF()

    # goes and reads file, sets up calc
    PSF.get_config()
    
    # now overwrite all config variables we want to change
    PSF.config.set_config()

    print(get_file_str(ind))




