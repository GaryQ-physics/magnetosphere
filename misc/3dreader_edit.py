
"""
Created on Sun Feb  2 09:58:51 2020

"""

import spacepy.pybats.bats as bt
import numpy as np
from matplotlib import pyplot as plt

# read in the 3d magnetosphere
filename = "3d__var_4_e20191012-131400-018.out"
data3d = bt.Bats2d(filename)

# look at keys:
print(data3d.keys())

# get the cell coordinates
x = data3d['x']
y = data3d['y']
z = data3d['z']
#bz = data3d['bz']

p, bx, by, bz =  np.array(data3d['p']), np.array(data3d['bx']), np.array(data3d['by']), np.array(data3d['bz'])

print(bx.size)
print(x.size)

plt.clf()

# plot p first:
ax1 = plt.subplot(121)
ax1.set_aspect('equal', 'box')
data3d.add_contour('x', 'y', 'bx', target = ax1)

plt.xlim([-50, 50])
plt.ylim([-50, 50])

plt.show()
