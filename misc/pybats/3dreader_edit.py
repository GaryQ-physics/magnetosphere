
"""
Created on Sun Feb  2 09:58:51 2020

Demonstrate reading the 3d .out datafiles using spacepy.pybats.bats, and then plotting with that data

Works for python 2 and 3
"""
import sys
import os
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../../' )
from config import conf

import spacepy.pybats.bats as bt
import numpy as np
from matplotlib import pyplot as plt
    
# read in the 3d magnetosphere
filename = conf['run_path'] + "3d__var_3_e20031120-070000-000.out"

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
