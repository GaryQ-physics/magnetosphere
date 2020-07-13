"""
pip install joblib

whern para is false:
    with N = 2560+1 (so dx=dy=dz=0.1) , it takes about 48s to run through each slice.
        Thus its estimated running to completion would take 34.1 hours.
        Assuming the 24 core comp in parallel decreases by factor 1/24, then 1.4hours.
    with N = 64+1 (so dx=dy=dz=4.), it prints Btot_EW = -28.53599693175405
    with N = 128+1 (so dx=dy=dz=2.), it prints Btot_EW = -13.94200518544455

Typical outputs:

    N = 65
    fullVolume = False
    slices = 33
    points in slice = 289
    total points = 9537
    X[0], X[-1], dx = -96.000000, 32.000000, 4.000000
    Y[0], Y[-1], dy = -32.000000, 32.000000, 4.000000
    Z[0], Z[-1], dz = -32.000000, 32.000000, 4.000000
    para = False
    time to process all slices (not including suming up) = 0.20150 min
    Btot = 
    [ -3.81502299   9.8473954  -40.06357733]
    Btot_norm = 41.43206276763158

    N = 129
    fullVolume = False
    slices = 65
    points in slice = 1089
    total points = 70785
    X[0], X[-1], dx = -96.000000, 32.000000, 2.000000
    Y[0], Y[-1], dy = -32.000000, 32.000000, 2.000000
    Z[0], Z[-1], dz = -32.000000, 32.000000, 2.000000
    para = True
    time to process all slices (not including suming up) = 0.15039 min
    Btot = 
    [-2.3255107  -1.41573583 -7.77186095]
    Btot_norm = 8.234933548023848

    N = 129
    fullVolume = True
    slices = 129
    points in slice = 16641
    total points = 2146689
    X[0], X[-1], dx = -224.000000, 32.000000, 2.000000
    Y[0], Y[-1], dy = -128.000000, 128.000000, 2.000000
    Z[0], Z[-1], dz = -128.000000, 128.000000, 2.000000
    para = True
    time to process all slices (not including suming up) = 0.42586 min
    Btot = 
    [-2.31410359 -1.34009783 -7.8129071 ]
    Btot_norm = 8.257872300049227

    N = 257
    fullVolume = False
    slices = 129
    points in slice = 4225
    total points = 545025
    X[0], X[-1], dx = -96.000000, 32.000000, 1.000000
    Y[0], Y[-1], dy = -32.000000, 32.000000, 1.000000
    Z[0], Z[-1], dz = -32.000000, 32.000000, 1.000000
    para = True
    time to process all slices (not including suming up) = 0.28954 min
    Btot = 
    [-4.93034021  0.05562719  0.85745892]
    Btot_norm = 5.00465630920011

    N = 257
    fullVolume = True
    slices = 257
    points in slice = 66049
    total points = 16974593
    X[0], X[-1], dx = -224.000000, 32.000000, 1.000000
    Y[0], Y[-1], dy = -128.000000, 128.000000, 1.000000
    Z[0], Z[-1], dz = -128.000000, 128.000000, 1.000000
    para = True
    time to process all slices (not including suming up) = 1.63295 min
    Btot = 
    [-4.92360583  0.13197522  0.81855364]
    Btot_norm = 4.992929186467655

    N = 2561
    fullVolume = False
    slices = 1281
    points in slice = 410881
    total points = 526338561
    X[0], X[-1], dx = -96.000000, 32.000000, 0.100000
    Y[0], Y[-1], dy = -32.000000, 32.000000, 0.100000
    Z[0], Z[-1], dz = -32.000000, 32.000000, 0.100000
    para = True
    time to process all slices (not including suming up) = 39.38515 min (26.29395 min with usekV)
    Btot = 
    [-3.13902583  3.25972791 -3.59094972]
    Btot_norm = 5.77704328474687

"""

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import numpy as np
import biot_savart as bs
import cxtransform as cx
from units_and_constants import phys
from probe import probe

############
Test = False
############


def run(time, mlat, mlon, para=True, fullVolume=False, n=2560, spacepy_like=False):
    N = n + 1

    if spacepy_like:
        D = 3.96875
        X = np.linspace(-D, D, 128)
        Y = np.linspace(-D, D, 128)
        Z = np.linspace(-D, D, 128)
    elif fullVolume:
        global_x_min = -224.
        global_x_max = 32.
        global_y_min = -128.
        global_y_max = 128.
        global_z_min = -128.
        global_z_max = 128.
        # difference of 256 for all

        X = np.linspace(global_x_min, global_x_max, N)
        Y = np.linspace(global_y_min, global_y_max, N)
        Z = np.linspace(global_z_min, global_z_max, N)
    else:
        global_x_min = -224.
        global_x_max = 32.
        global_y_min = -128.
        global_y_max = 128.
        global_z_min = -128.
        global_z_max = 128.
        # difference of 256 for all

        X = np.linspace(global_x_min + 128., global_x_max, n/2+1)
        Y = np.linspace(global_y_min/4, global_y_max/4, n/4+1)
        Z = np.linspace(global_z_min/4, global_z_max/4, n/4+1)


    dx = X[1]-X[0]
    dy = Y[1]-Y[0]
    dz = Z[1]-Z[0]

    if spacepy_like:
        assert(dx == 0.0625)


    Gy, Gz = np.meshgrid(Y,Z)
    Gy = Gy.flatten(order='C')
    Gz = Gz.flatten(order='C')

    x0 = cx.MAGtoGSM([1., mlat, mlon], time, 'sph', 'car')
    print(x0)
    #Npole = cx.GEOtoGSM([0., 0., 1.], time, 'car', 'car')

    def dBslice(i, debug=False):
        Grid = np.column_stack([X[i]*np.ones(Gy.shape), Gy, Gz])
        J_kameleon = probe(time, Grid, ['jx','jy','jz'])
        J = J_kameleon*(phys['muA']/phys['m']**2)
        if debug:
            print(Grid.shape)
            print(J.shape)

        deltaB = bs.deltaB('deltaB', x0, Grid, J, V_char = dx*dy*dz)
        return deltaB
        #return bs.B_EW(X0, Grid, J, Npole, dx*dy*dz)

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
        B_slices = Parallel(n_jobs=num_cores)(delayed(dBslice)(j) for j in range(X.size))
        B_slices = np.array(B_slices) # was list of (3,) numpy arrays
    else:
        B_slices = np.nan*np.empty((X.size, 3))
        for j in range(X.size):
            B_slices[j,:] = dBslice(j)

    tf = t_module.time()
    #print(B_slices.shape)

    print('N = ' + str(N))
    print('fullVolume = ' + str(fullVolume))
    print('slices = ' + str(X.size))
    print('points in slice = ' + str(Y.size*Z.size))
    print('total points = ' + str(X.size*Y.size*Z.size))
    print('X[0], X[-1], dx = {0:f}, {1:f}, {2:f}'.format(X[0], X[-1], dx))
    print('Y[0], Y[-1], dy = {0:f}, {1:f}, {2:f}'.format(Y[0], Y[-1], dy))
    print('Z[0], Z[-1], dz = {0:f}, {1:f}, {2:f}'.format(Z[0], Z[-1], dz))
    print('para = ' + str(para))
    print('time to process all slices (not including suming up) = {0:.5f} min'\
            .format((tf-to)/60.))

    Btot = np.sum(B_slices, axis=0)

    print('Btot = \n' + str(Btot))
    print('Btot_norm = ' + str(np.linalg.norm(Btot)))

    return Btot


if Test:
    data = np.array([[2003, 11, 20, 7, 0, 57.50, 176.00]])
    #data = events()

    time = data[0, 0:5]
    mlat = data[0, 5]
    mlon = data[0, 6]
    Event = data[0, :]

    para = True
    fullVolume = False

    n = 128 # for testing
    #n = 2560
    N = n + 1 # for testing

    run(time, mlat, mlon, para=para, fullVolume=fullVolume, n=n)
