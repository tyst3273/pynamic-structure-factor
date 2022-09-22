
import matplotlib.pyplot as plt
import numpy as np

_ = np.loadtxt('sqw')
_en = _[0,:]
_n = _[1,:]
_eo, _o = np.loadtxt('tmp/old_sqw',unpack=True)

print(_en.shape)
print(_n.shape)
print(_eo.shape)
print(_o.shape)

plt.plot(_en,_n,c='r')
plt.plot(_eo,_o,c='b')
plt.show()

