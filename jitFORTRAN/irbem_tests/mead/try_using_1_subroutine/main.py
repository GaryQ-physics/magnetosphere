import numpy as np
import pymod

try:
    print(pymod.meadv)
    print('yes pymod.meadv')
except:
    print('no pymod.meadv')

try:
    print(pymod.mead)
    print('yes pymod.mead')
except:
    print('no pymod.mead')


N = 10
x = np.arange(10, dtype=np.float64)+1.
y = np.arange(10, dtype=np.float64)+1.
z = np.arange(10, dtype=np.float64)+1.

bx = np.nan*np.empty((10,), dtype=np.float64)
by = np.nan*np.empty((10,), dtype=np.float64)
bz = np.nan*np.empty((10,), dtype=np.float64)

pymod.meadv(x, y, z, 1, bx, by, bz, N)

print(bx, by, bz)
