import sys
import os
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../../' )
from config import conf
import numpy as np


import spacepy.pybats.bats as bats

from util import time2filename
data = np.array([[2003, 11, 20, 7, 0, 57.50, 176.00]])
#data = events()
time = data[0, 0:5]
mlat = data[0, 5]
mlon = data[0, 6]
Event = data[0, :]

# read in the 3d magnetosphere
filename = time2filename(time, extention='.out')
data3d = bats.Bats2d(filename)

from units_and_constants import phys

# get the cell coordinates
x = data3d['x']
y = data3d['y']
z = data3d['z']
jx = data3d['jx']
jy = data3d['jy']
jz = data3d['jz']

x = np.array(x, dtype=float)
y = np.array(y, dtype=float)
z = np.array(z, dtype=float)
jx = np.array(jx, dtype=float)
jy = np.array(jy, dtype=float)
jz = np.array(jz, dtype=float)

D = 3.96875
xax = np.linspace(-D, D, 128)
yax = np.linspace(-D, D, 128)
zax = np.linspace(-D, D, 128)
dx = 0.0625
dV = dx**3

assert((D+D)/(128-1) == dx)

Int = 5
Slice = x==(0.03125+Int*0.0625)

Tr = np.all([-D<=x, x<=D, -D<=y, y<=D, -D<=z, z<=D, Slice], axis=0)
x_ = x[Tr]
y_ = y[Tr]
z_ = z[Tr]
jx_ = jx[Tr]
jy_ = jy[Tr]
jz_ = jz[Tr]

J_un = np.column_stack([jx_, jy_, jz_])
X = np.column_stack([x_, y_, z_])

from probe import probe
J_kameleon = probe(time, X, var=['jx', 'jy', 'jz'], usekV=False)

print(J_un.shape)
print(J_kameleon.shape)

"""
you can see in the resulting graph that the kameleon exactly agrees with
spacepy except on a few right at the boundary of this  -3.96875, 3.96875
cube region
"""
tru = np.all([jx_ == J_kameleon[:,0],\
              jy_ == J_kameleon[:,1],\
              jz_ == J_kameleon[:,2]], axis=0)
fal = np.logical_not(tru)


y_true = y_[tru]
z_true = z_[tru]
y_false = y_[fal]
z_false = z_[fal]

import matplotlib.pyplot as plt
plt.plot(y_true, z_true, marker='.', color='b', linestyle='')
plt.plot(y_false, z_false, marker='.', color='r', linestyle='', markersize=2.)
plt.xlabel('xlable')
plt.ylabel('ylable')
plt.show()
