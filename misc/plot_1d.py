# plot_1d

import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import _CCMC as ccmc
from units_and_constants import phys

from probe import probe, probe_vect


time = (2003, 11, 20, 7, 0, 0)
ds = 0.05
u = np.array([1., 0., 0.])

s = np.arange(0., 3. + ds, ds)

U = np.repeat([u], s.size, axis=0)

line = U*s[:,np.newaxis]

print(line.shape)
print(line)


Jin = probe_vect(time, line, 'j')*(phys['muA']/phys['m']**2)

print(Jin.shape)

import matplotlib.pyplot as plt
plt.plot(s, Jin[:,2], 'r')
plt.xlabel('s / R_e (length along line in direction ({0:.2f}, {1:.2f}, {2:.2f}))'.format(*tuple(u)))
plt.ylabel('jy / (muA/R_e**2)')
plt.axvline(x=1.25)
plt.axvline(x=1.2)
plt.show()
