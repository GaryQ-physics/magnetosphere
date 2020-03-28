import spacepy.pybats.bats as bats
import spacepy.pybats.bats as bt
import numpy as np
from matplotlib import pyplot as plt
############################################################################
# read in the 3d magnetosphere
filename = "3d__var_4_e20191012-131400-018.out"
data3d = bats.Bats2d(filename)

# look at keys:
print(data3d.keys())

# get the cell coordinates
x = data3d['x']
y = data3d['y']
z = data3d['z']
bz = data3d['bz']










