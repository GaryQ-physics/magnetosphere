# plot_1d

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
from units_and_constants import phys
from probe import probe, probe_vect

time = (2003, 11, 20, 7, 0, 0)
ds = 0.05
u = np.array([1., 0., 0.])

s = np.arange(0., 3. + ds, ds)

U = np.repeat([u], s.size, axis=0)

line = U*s[:, np.newaxis]

print(line.shape)
print(line)

Jin = probe_vect(time, line, 'j')*(phys['muA']/phys['m']**2)

print(Jin.shape)

plt.plot(s, Jin[:,2], 'r')
plt.xlabel('$s/R_E$ (length along line in direction ({0:.2f}, {1:.2f}, {2:.2f}))'.format(*tuple(u)))
plt.ylabel('$J_y/(\mu A/R_E^2)$')
plt.axvline(x=1.25)
plt.axvline(x=1.2)
plt.grid()
plt.show()
