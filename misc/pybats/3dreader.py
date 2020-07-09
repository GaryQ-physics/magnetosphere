"""
Demonstrate reading the 3d .out datafiles using spacepy.pybats.bats

Works for python 2 and 3
"""

import sys
import os
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../../' )
from config import conf

import spacepy.pybats.bats as bats

############################################################################
# read in the 3d magnetosphere
filename = conf['run_path'] + "3d__var_3_e20031120-070000-000.out"
data3d = bats.Bats2d(filename)

# look at keys:
print(data3d.keys())

import numpy as np
# get the cell coordinates
x = data3d['x']
y = data3d['y']
z = data3d['z']
bz = data3d['bz']
jy = data3d['jy']
jz = data3d['jz']
p = data3d['p']

x = np.array(x, dtype=float)
y = np.array(y)
z = np.array(z)
p = np.array(p)

#print(np.sqrt(np.min(x**2 + y**2 + z**2)))
Int = 58 # 62
Tr = x == (0.03125 + Int*0.0625)

eps = 0.011
a = 0.04125
#Tr = np.logical_and(a - eps <= x, x <= a + eps)
#Tr = np.logical_and(y==0.03125, x==0.03125)
x_ = x[Tr]
y_ = y[Tr]
z_ = z[Tr]
p_ = p[Tr]
jy_ = jy[Tr]
jz_ = jz[Tr]

#print(x_)
#print(y_)
#print(z_)
print(np.min(z_), np.max(z_))
print(z_.shape)

yax = np.linspace(-3.96875, 3.96875, 128)
zax = np.linspace(-3.96875, 3.96875, 128)

#print(yax)

Y,Z = np.meshgrid(yax, zax)

import matplotlib.pyplot as plt
plt.plot(y_, z_, marker='.', color='b', linestyle='')
plt.plot(Y.flatten(), Z.flatten(), marker='.', color='r', linestyle='', markersize=2.)
plt.xlabel('xlable')
plt.ylabel('ylable')
#plt.axvline(x=1.25)
#plt.axvline(x=1.2)
plt.show()
