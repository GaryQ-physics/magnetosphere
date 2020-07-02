import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import spacepy.pybats.bats as bats
#import _CCMC as ccmc

from units_and_constants import phys
from probe import probe

time = (2003, 11, 20, 7, 0, 0)
# read in the 3d magnetosphere
filename = conf['run_path'] + "3d__var_3_e20031120-070000-000.out"
data3d = bats.Bats2d(filename)

# get the cell coordinates
x = data3d['x']
y = data3d['y']
z = data3d['z']
bz = data3d['bz']
jy = data3d['jy']
jz = data3d['jz']
p = data3d['p']

#Tr = x==0.03125
x = np.array(x, dtype=float)
y = np.array(y, dtype=float)
z = np.array(z, dtype=float)


eps = 0.011
a = 0.04125
Tr = np.logical_and(a - eps <= x, x <= a + eps)
Tr = np.logical_and(y==0.03125, z==0.03125)
x_ = x[Tr]
y_ = y[Tr]
z_ = z[Tr]
p_ = p[Tr]
jy_ = jy[Tr]
jz_ = jz[Tr]




line = np.column_stack([x_, y_, z_])


Jin = probe(time, line, ['jx','jy','jz'])*(phys['muA']/phys['m']**2)



import matplotlib.pyplot as plt
plt.plot(x_, jz_, '.')
#plt.plot(line[:,0], Jin[:,2], '.')
plt.xlabel('xlable')
plt.ylabel('ylable')
plt.axvline(x=1.25)
plt.axvline(x=1.2)
plt.show()
