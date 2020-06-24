"""
with N = 2560+1 (so dx=dy=dz=0.1) , it takes about 48s to run through each slice.
    Thus its estimated running to completion would take 34.1 hours.
    Assuming the 24 core comp in parallel decreases by factor 1/24, then 1.4hours.

with N = 64+1 (so dx=dy=dz=4.), it prints Btot_EW = -28.53599693175405

with N = 128+1 (so dx=dy=dz=2.), it prints Btot_EW = -13.94200518544455

"""

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import numpy as np
import biot_savart as bs
import pos_sun as ps
from units_and_constants import phys
#from util import time2filename, filemeta
from probe import probe_vect


line_list = [2003, 11, 20, 7, 0, 176.00, 57.50]
time = line_list[0:5] + [0, 0.]
Event = time + line_list[5:7]
MLON = Event[5]
MLAT = Event[6]

Test = False

rbody = 1.25
global_x_min = -224.
global_x_max = 32.
global_y_min = -128.
global_y_max = 128.
global_z_min = -128.
global_z_max = 128.
# difference of 256 for all
diff = 256.

N = 128 + 1 # for testing
#N = 2560 + 1 
dx = diff/(N-1)
dy = diff/(N-1)
dz = diff/(N-1)

print('slices =', N)
print('points in slice =', N**2)
print('total points =', N**3)

X = np.linspace(global_x_min, global_x_max, N)
Y = np.linspace(global_y_min, global_y_max, N)
Z = np.linspace(global_z_min, global_z_max, N)

Gy, Gz = np.meshgrid(Y,Z)
Gy = Gy.flatten(order='C')
Gz = Gz.flatten(order='C')

if Test:
    X0 = np.array([0., 0.75, 0.])
    Npole = np.array([0., 0., 1.])
else:
    X0 = ps.MAGtoGSM([1., MLAT, MLON], time[0:6], 'sph', 'car')
    Npole = ps.GEOtoGSM([0., 0., 1.], time[0:6], 'car', 'car')

def dBslice(i, debug=False):
    Grid = np.column_stack([X[i]*np.ones(Gy.shape), Gy, Gz])
    J_kameleon = probe_vect(time, Grid, 'j')
    J = J_kameleon*(phys['muA']/phys['m']**2)
    if debug:
        print(Grid.shape)
        print(J.shape)

    return bs.B_EW(X0, Grid, J, Npole, dx*dy*dz)

dB_EW = np.nan*np.empty(X.shape)
for j in range(X.size):
    dB_EW[j] = dBslice(j)

Btot_EW = np.sum(dB_EW)

print(dB_EW.shape)
print(Btot_EW)

