"""
pip install joblib

whern para is false:
    with N = 2560+1 (so dx=dy=dz=0.1) , it takes about 48s to run through each slice.
        Thus its estimated running to completion would take 34.1 hours.
        Assuming the 24 core comp in parallel decreases by factor 1/24, then 1.4hours.
    with N = 64+1 (so dx=dy=dz=4.), it prints Btot_EW = -28.53599693175405
    with N = 128+1 (so dx=dy=dz=2.), it prints Btot_EW = -13.94200518544455

Typical outputs:

    Nx, Ny, Nz = 65, 33, 33
    fullVolume = False
    spacepy_like = False
    slices = 65
    points in slice = 1089
    total points = 70785
    X[0], X[-1], dx = -96.000000, 32.000000, 2.000000
    Y[0], Y[-1], dy = -32.000000, 32.000000, 2.000000
    Z[0], Z[-1], dz = -32.000000, 32.000000, 2.000000
    para = True
    time to process all slices (not including suming up) = 0.13764 min
    Btot = 
    [-2.3255107  -1.41573583 -7.77186095]
    Btot_norm = 8.234933548023848

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


def run(time, mlat, mlon, para=True,
        fullVolume=False, spacepy_like=False, 
        N=None, L=None, eps=None,
        Nx=None, xlims=None, dx=None,
        Ny=None, ylims=None, dy=None,
        Nz=None, zlims=None, dz=None,
        print_output=False, tolerance=1e-13):

###### make X, Y, and Z ##############
    assert(not (fullVolume and spacepy))

    if spacepy_like:
        L = 3.96875
        eps = 0.0625
        N = 128
        assert((L+L)/(N-1) == eps)
        assert((2*L)/(N-1) == eps)
        tolerance = 0.
    if fullVolume:
        xlims = (-224., 32.)
        ylims = (-128., 128.)
        zlims = (-128., 128.)

    if L != None:
        xlims = (-L, L)
        ylims = (-L, L)
        zlims = (-L, L)

    if N != None:
        Nx = N
        Ny = N
        Nz = N

    if eps != None:
        dx = eps
        dy = eps
        dz = eps


    assert(xlims!=None and xlims!=None and zlims!=None)

    if Nx != None:
        X = np.linspace(xlims[0], xlims[1], Nx)
    elif dx != None:
        X = np.arange(xlims[0], xlims[1] + dx, dx)
    else:
        assert(False)

    if Ny != None:
        Y = np.linspace(ylims[0], ylims[1], Ny)
    elif dy != None:
        Y = np.arange(ylims[0], ylims[1] + dy, dy)
    else:
        assert(False)

    if Nz != None:
        Z = np.linspace(zlims[0], zlims[1], Nz)
    elif dz != None:
        Z = np.arange(zlims[0], zlims[1] + dz, dz)
    else:
        assert(False)

    dx_check = X[1]-X[0]
    dy_check = Y[1]-Y[0]
    dz_check = Z[1]-Z[0]
    x_range_check = X[-1]-X[0]
    y_range_check = Y[-1]-Y[0]
    z_range_check = Z[-1]-Z[0]
    Nx_check = X.size
    Ny_check = Y.size
    Nz_check = Z.size

    assert(Nx_check == Nx or Nx == None)
    if dx != None:
        assert(np.abs(dx_check - dx) <= tolerance)
    assert(np.abs(xlims[1] - xlims[0] - x_range_check) <= tolerance)

    assert(Ny_check == Ny or Ny == None)
    if dy != None:
        assert(np.abs(dy_check - dy) <= tolerance or dy == None)
    assert(np.abs(ylims[1] - ylims[0] - y_range_check) <= tolerance)

    assert(Nz_check == Nz or Nz == None)
    if dz != None:
        assert(np.abs(dz_check - dz) <= tolerance or dz == None)
    assert(np.abs(zlims[1] - zlims[0] - z_range_check) <= tolerance)

######################################

    dx = dx_check
    dy = dy_check
    dz = dz_check

    Nx = Nx_check
    Ny = Ny_check
    Nz = Nz_check


    Gy, Gz = np.meshgrid(Y,Z)
    Gy = Gy.flatten(order='C')
    Gz = Gz.flatten(order='C')

    x0 = cx.MAGtoGSM([1., mlat, mlon], time, 'sph', 'car')
    if print_output:
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

    if print_output:
        print('Nx, Ny, Nz = {0:d}, {1:d}, {2:d}'.format(Nx,Ny,Nz))
        print('fullVolume = ' + str(fullVolume))
        print('spacepy_like = ' + str(spacepy_like))
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
    if print_output:
        print('Btot = \n' + str(Btot))
        print('Btot_norm = ' + str(np.linalg.norm(Btot)))

    return Btot


if Test:
    global_x_min = -224.
    global_x_max = 32.
    global_y_min = -128.
    global_y_max = 128.
    global_z_min = -128.
    global_z_max = 128.

    data = np.array([[2003, 11, 20, 7, 0, 57.50, 176.00]])
    #data = events()

    time = data[0, 0:5]
    mlat = data[0, 5]
    mlon = data[0, 6]
    Event = data[0, :]

    para = True
    fullVolume = False
    spacepy_like = False


    n = 128 # for testing
    #n = 2560

    run(time, mlat, mlon, para=para,
        fullVolume=fullVolume, spacepy_like=spacepy_like,
        Nx=n/2+1, xlims=(global_x_min + 128., global_x_max), dx=None,
        Ny=n/4+1, ylims=(global_y_min/4, global_y_max/4), dy=None,
        Nz=n/4+1, zlims=(global_z_min/4, global_z_max/4), dz=None,
        print_output=True)
