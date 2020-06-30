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

from probe import probe

#RUN PARAMETERS
debug = False
time = (2003, 11, 20, 7, 0, 0)
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

Jin = probe(time, line, ['jx','jy','jz'])*(phys['muA']/phys['m']**2)

if debug:
    print(line.shape)
    print(line)
    print(Jin.shape)

plt.plot(s, Jin[:,2], 'r')
plt.xlabel('$s / R_e$ (length along direction \
    ({0:.2f}, {1:.2f}, {2:.2f}))'.format(*tuple(u)))
plt.ylabel('$J_y / (\mu A / R_e^2)$')
plt.axvline(x=1.25)
plt.axvline(x=1.2)
plt.grid()
plt.show()
