"""
Creates a plot of a kameleon variable along a line trough the origing
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
from probe import GetRunData
import spacepy.pybats.bats as bats
import util


#RUN PARAMETERS
debug = False
run = 'DIPTSUR2'
time = (2019,9,2,6,30,0,0)

D = 3.96875
delta = 0.0625

# read in the 3d magnetosphere
data3d = bats.Bats2d(util.time2CDFfilename(run, time)[:-4])

# get outfile data
x = data3d['x']
y = data3d['y']
z = data3d['z']
jx = data3d['jx']
jy = data3d['jy']
jz = data3d['jz']

x = np.array(x, dtype=float)
y = np.array(y, dtype=float)
z = np.array(z, dtype=float)

# restrict to line
Tr = np.logical_and(y == (delta/2.), z == (delta/2.))
x_Tr = x[Tr]
y_Tr = y[Tr]
z_Tr = z[Tr]
jz_Tr = jz[Tr]

# make finer line
n = 10

xax = np.linspace(-D, D, n*128)
line = np.column_stack([xax, (delta/2.)*np.ones(n*128), (delta/2.)*np.ones(n*128)])

#camelon interpolate
Jz = GetRunData(run, time, line, 'jz')



plt.plot(x_Tr , jz_Tr, marker='x', color='b', linestyle='')
plt.plot(xax, Jz, marker='+', color='r', linestyle='',)
#plt.xlabel('$s / R_e$ (length along direction ({0:.2f}, {1:.2f}, {2:.2f}))'.format(*tuple()))
#plt.ylabel('$J_z$')
#plt.axvline(x=1.25)
#plt.axvline(x=1.2)
plt.grid()
plt.show()
