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
x = np.array(x)
y = np.array(y)
z = np.array(z)

print(np.sqrt(np.min(x**2 + y**2 + z**2)))
