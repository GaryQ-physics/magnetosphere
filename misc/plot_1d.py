"""
Creates a plot plot of a kameleon variable along a line trough the origing
    on the x-axis is (signed) distance from the origin of a point along the line
    on the y axis is the kameleon variable at that point
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
from units_and_constants import phys


#RUN PARAMETERS
debug = False
time = (2003, 11, 20, 7, 0, 0)
filename = conf['run_path'] + "3d__var_3_e20031120-070000-000.out"

from probe import probe
import spacepy.pybats.bats as bats

# read in the 3d magnetosphere
data3d = bats.Bats2d(filename)

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


D = 3.96875
delta = 0.0625
xax = np.linspace(-D, D, 128)
yax = np.linspace(-D, D, 128)
zax = np.linspace(-D, D, 128)

if False:
    # direction vector of line
    u = np.array([1., 0., 0.]) # along x-axis
    # max length from origin that is used
    L = 3.
    # step size of length along line
    ds = 0.05

    s = np.arange(0., L + ds, ds)
    # Nx3 array where each row is a duplicate of u
    U = np.repeat([u], s.size, axis=0) 
    # Nx3 array where i'th row is i'th row vector of U (so just u)
    #     multiplied by i'th element of s
    line = U*s[:,np.newaxis]

line = np.column_stack([xax, (delta/2.)*np.ones(128), (delta/2.)*np.ones(128)])

Tr = np.logical_and(y == (delta/2.), z == (delta/2.))
x_ = x[Tr]
y_ = y[Tr]
z_ = z[Tr]
claim = np.column_stack([x_, y_, z_]) == line
jz_ = jz[Tr]

assert(np.all(claim))

J_kam = probe(time, line, ['jx','jy','jz'])

if debug:
    print(claim)
    print(claim.shape)
    print(line.shape)
    print(line)
    print(x_)
    print(Jin.shape)

plt.plot(x_ , jz_, marker='.', color='b', linestyle='')
plt.plot(line[:,0], J_kam[:,2], marker='.', color='r', linestyle='', markersize=2.)
#plt.xlabel('$s / R_e$ (length along direction ({0:.2f}, {1:.2f}, {2:.2f}))'.format(*tuple()))
plt.ylabel('$J_z$')
plt.axvline(x=1.25)
plt.axvline(x=1.2)
plt.grid()
plt.show()
