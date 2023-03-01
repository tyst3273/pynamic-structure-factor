
from PSF import c_PSF

input_file = 'input_cmap.py'
PSF = c_PSF(input_file)

print(PSF.__dir__())

PSF.setup_calculation(num_trajectory_blocks=10)
PSF.run()

