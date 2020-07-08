import os
import sys
import numpy as np

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../../' )
from config import conf

"""
Compare Kameleon-computed interpolation values at grid points derived from
original data file (read using PyBats module in SpacePy). Differences are
zero except at outer boundary of high-resolution grid where the differences
are ~10000 times smaller than np.finfo(np.float32).eps and ~100000 times
larger than np.finfo(np.float64).eps. (The PyBats file reader returns data 
of type float32). The differences are likely explained by casting
(the Kamelon library takes input of float64 and returns float64). However,
the fact that differences occur only on the outer edge is suspicious as is
the fact that the differences are not on the order of np.finfo(np.float32).eps.

Note:
    np.finfo(np.float32).eps == 1.1920929e-07

    np.finfo(np.float64).eps == 2.220446049250313e-16
    np.finfo(float).eps      == 2.220446049250313e-16
"""

import spacepy.pybats.bats as bats

from util import time2filename
time = [2003, 11, 20, 7, 0]

# read original data file
filename = time2filename(time, extension='.out')
data3d = bats.Bats2d(filename)

# get the cell coordinates and associated values of J components
x = data3d['x']
y = data3d['y']
z = data3d['z']
jx = data3d['jx']
jy = data3d['jy']
jz = data3d['jz']

x = np.array(x, dtype=float)  #passing to probe() does not work unless specify dtype
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
x_cut = .03125+Int*0.0625
Slice = x == x_cut

Tr = np.all([-D<=x, x<=D, -D<=y, y<=D, -D<=z, z<=D, Slice], axis=0)
x_ = x[Tr]
y_ = y[Tr]
z_ = z[Tr]
jx_ = jx[Tr]
jy_ = jy[Tr]
jz_ = jz[Tr]

X = np.column_stack([x_, y_, z_])

from probe import probe
J_kameleon = probe(time, X, var=['jx', 'jy', 'jz'], usekV=False)

print(jx_.shape)
print(J_kameleon.shape)

tru = np.all([jx_ == J_kameleon[:,0],
              jy_ == J_kameleon[:,1],
              jz_ == J_kameleon[:,2]], axis=0)
fal = np.logical_not(tru)


y_true = y_[tru]
z_true = z_[tru]
y_false = y_[fal]
z_false = z_[fal]

import matplotlib.pyplot as plt
plt.figure(dpi=96*3)
plt.title('X = ' + str(x_cut))
plt.plot(y_true, z_true, marker='.', color='b', linestyle='', markersize=0.05)
plt.plot(y_false, z_false, marker='.', color='r', linestyle='', markersize=1)
plt.axis('square')
plt.xlabel('Y')
plt.ylabel('Z')
plt.show()
plt.savefig('check_kameleon1.png', dpi=96*3)

plt.figure(dpi=96*3)
plt.plot((jx_ - J_kameleon[:,0])/np.finfo(np.float32).eps)
plt.xlabel('index')
plt.title('($j_x$ actual - $j_x$ kameleon)/np.finfo(np.float32).eps')
plt.show()
plt.savefig('check_kameleon2.png', dpi=96*3)

