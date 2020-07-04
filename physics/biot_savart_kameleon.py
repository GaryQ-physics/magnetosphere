"""
pip install joblib

whern para is false:
    with N = 2560+1 (so dx=dy=dz=0.1) , it takes about 48s to run through each slice.
        Thus its estimated running to completion would take 34.1 hours.
        Assuming the 24 core comp in parallel decreases by factor 1/24, then 1.4hours.
    with N = 64+1 (so dx=dy=dz=4.), it prints Btot_EW = -28.53599693175405
    with N = 128+1 (so dx=dy=dz=2.), it prints Btot_EW = -13.94200518544455

Typical outputs:
    N = 129
    fullVolume = True
    slices = 129
    points in slice = 16641
    total points = 2146689
    X[0], X[-1], dx = -224.000000, 32.000000, 2.000000
    Y[0], Y[-1], dy = -128.000000, 128.000000, 2.000000
    Z[0], Z[-1], dz = -128.000000, 128.000000, 2.000000
    para = True
    time to process all slices (not including suming up) = 23.13155 sec
    Btot_EW = -13.94200518544455

    N = 129
    fullVolume = False
    slices = 65
    points in slice = 1089
    total points = 70785
    X[0], X[-1], dx = -96.000000, 32.000000, 2.000000
    Y[0], Y[-1], dy = -32.000000, 32.000000, 2.000000
    Z[0], Z[-1], dz = -32.000000, 32.000000, 2.000000
    para = True
    time to process all slices (not including suming up) = 8.25513 sec
    Btot_EW = -13.862501156128777

    N = 2561
    fullVolume = False
    slices = 1281
    points in slice = 410881
    total points = 526338561
    X[0], X[-1], dx = -96.000000, 32.000000, 0.100000
    Y[0], Y[-1], dy = -32.000000, 32.000000, 0.100000
    Z[0], Z[-1], dz = -32.000000, 32.000000, 0.100000
    para = True
    time to process all slices (not including suming up) = 2348.29029 sec
    Btot_EW = -8.14956342968862

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
from probe import probe


line_list = [2003, 11, 20, 7, 0, 176.00, 57.50]
time = line_list[0:5] + [0, 0.]
Event = time + line_list[5:7]
MLON = Event[5]
MLAT = Event[6]

Test = False
para = True
fullVolume = False

rbody = 1.25
global_x_min = -224.
global_x_max = 32.
global_y_min = -128.
global_y_max = 128.
global_z_min = -128.
global_z_max = 128.
# difference of 256 for all
diff = 256.

#n = 128 # for testing
n = 2560
N = n + 1 # for testing


#dx = diff/(N-1)
#dy = diff/(N-1)
#dz = diff/(N-1)

if fullVolume:
    X = np.linspace(global_x_min, global_x_max, N)
    Y = np.linspace(global_y_min, global_y_max, N)
    Z = np.linspace(global_z_min, global_z_max, N)
else:
    X = np.linspace(global_x_min + 128., global_x_max, n/2+1)
    Y = np.linspace(global_y_min/4, global_y_max/4, n/4+1)
    Z = np.linspace(global_z_min/4, global_z_max/4, n/4+1)

dx = X[1]-X[0]
dy = Y[1]-Y[0]
dz = Z[1]-Z[0]

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
    J_kameleon = probe(time, Grid, ['jx','jy','jz'])
    J = J_kameleon*(phys['muA']/phys['m']**2)
    if debug:
        print(Grid.shape)
        print(J.shape)

    return bs.B_EW(X0, Grid, J, Npole, dx*dy*dz)

'''
if para:
    from joblib import Parallel, delayed
    import multiprocessing
    num_cores = multiprocessing.cpu_count()
    if num_cores is not None and num_cores > len(vars):
        num_cores = len(vars)
    print('Parallel processing {0:d} variable(s) using {1:d} cores'\
          .format(len(opts.keys()), num_cores))
    results = Parallel(n_jobs=num_cores)(\
                delayed(process_var)(var, opts[var]) for var in vars)
else:
    print('Serial processing {0:d} variable(s).'.format(len(opts.keys())))
    i = 0
    results = []
    for var in vars:
        results.append(process_var(var, opts[var]))
        i = i + 1

return results
'''

import time as t_module
to = t_module.time()

if para:
    from joblib import Parallel, delayed
    import multiprocessing
    num_cores = multiprocessing.cpu_count()
    if num_cores is not None and num_cores > X.size:
        num_cores = X.size
    print('Parallel processing {0:d} slices(s) using {1:d} cores'\
          .format(X.size, num_cores))
    dB_EW = Parallel(n_jobs=num_cores)(delayed(dBslice)(j) for j in range(X.size))
    dB_EW = np.array(dB_EW) # was list
else:
    dB_EW = np.nan*np.empty(X.shape)
    for j in range(X.size):
        dB_EW[j] = dBslice(j)

tf = t_module.time()

print('N = ' + str(N))
print('fullVolume = ' + str(fullVolume))
print('slices = ' + str(X.size))
print('points in slice = ' + str(Y.size*Z.size))
print('total points = ' + str(X.size*Y.size*Z.size))
print('X[0], X[-1], dx = {0:f}, {1:f}, {2:f}'.format(X[0], X[-1], dx))
print('Y[0], Y[-1], dy = {0:f}, {1:f}, {2:f}'.format(Y[0], Y[-1], dy))
print('Z[0], Z[-1], dz = {0:f}, {1:f}, {2:f}'.format(Z[0], Z[-1], dz))
print('para = ' + str(para))
print('time to process all slices (not including suming up) = {0:.5f} sec'\
        .format(tf-to))

Btot_EW = np.sum(dB_EW)

#print(dB_EW.shape)
print('Btot_EW = ' + str(Btot_EW))
