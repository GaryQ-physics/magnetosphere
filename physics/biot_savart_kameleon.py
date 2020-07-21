"""
pip install joblib

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
from util import tpad, time2filename


def run(time, mlat, mlon, filename=None, para=True,
        fullVolume=False, fineVolume=False,
        xlims=(-56., 8.), ylims=(-32., 32.), zlims=(-32., 32.),
        d=0.125,
        dx=None, dy=None, dz=None,
        N=None,
        Nx=None, Ny=None, Nz=None,
        L=None,
        print_output=False, tolerance=1e-13):
    """

    spacepy_like -> native fine grid spacing

    Returns (Btot)
    -------
    (3,) numpy array
    
        using kameleon, through probe, returns the result of biot savart 
        integration on regular grid given by a mesh of 1d arrays X,Y,Z 
        (specified by input parameters) for the field at x0 corresponding to 
        mlat, mlon on earth.

    
    Parameters
    ----------
    time : tuple/list/array
        the time at which probe will evaluate the current.
    mlat : float
        magentic latitude of x0 (radius is 1 R_e).
    mlon : float
        magentic latitude of x0 
    para : boolean, OPTIONAL
        if true, then computation run in parallel on multiple processors. Output
        returned is identical.
        The default is True.
    Nx : integer, optional
        the number of points in X 
        The default is None.
    xlims : list/tuple, optional
        list/tuple of length 2
        X ranges from X[0]=xlims[0] to X[-1]=xlims[1]
        The default is None.
    dx : float, optional
        stepsize between points in X 
        The default is None.
    Ny : TYPE, optional
        DESCRIPTION. The default is None.
    ylims : TYPE, optional
        DESCRIPTION. The default is None.
    dy : TYPE, optional
        DESCRIPTION. The default is None.
    Nz : TYPE, optional
        DESCRIPTION. The default is None.
    zlims : TYPE, optional
        DESCRIPTION. The default is None.
    dz : TYPE, optional
        DESCRIPTION. The default is None.
    N : integer, optional
        if not None, it overrides Nx, Ny, and Nz to all be N
        The default is None.
    L : float, optional
        if not None, it overrides xlims,ylims,and zlims to all be (-L,L)
        The default is None.
    d : float, optional
        if not None, it overrides dx, dy, and dz to all be d
        The default is None.
    fullVolume : boolean, OPTIONAL
        if True, it overides xlims, ylims and zlims to be the limits of the 
        full kameleon grid. 
        The default is False.
    spacepy_like : bolean, OPTIONAL
        if True, it overides everything so that the grid is like the fine mesh 
        grid in the SWMF files except extend to be regular over a whole cube. 
        The default is False.
    print_output : boolean, optional
        prints output on command line. 
        The default is False.
    tolerance : floar, optional
        the tolerance to which which N and d type choices must be consistent.
        slight errors can be due to floating point arithmetic. 
        The default is 1e-13.
        Note if spacepy_like=True, this is overidden to 0.

    """
    ###### make X, Y, and Z ###########################
    assert(not (fullVolume and spacepy))

    if fineVolume:
        L = 3.96875
        d = 0.0625
        N = 128
        assert((L+L)/(N-1) == d)
        assert((2*L)/(N-1) == d)
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

    if d != None:
        dx = d
        dy = d
        dz = d


    assert(xlims!=None and ylims!=None and zlims!=None)

    if Nx != None:
        X = np.linspace(xlims[0], xlims[1], Nx)
    elif dx != None:
        X = np.arange(xlims[0], xlims[1] + dx, dx)
        #X = np.arange(xlims[0], xlims[1], dx, endpoint=True) doesnt work
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

    ############################################

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

    if filename == None:
        filename = time2filename(time)

    def dBslice(i, debug=False):
        Grid = np.column_stack([X[i]*np.ones(Gy.shape), Gy, Gz])
        J_kameleon = probe(filename, Grid, ['jx','jy','jz'], usekV=True)
        J = J_kameleon*(phys['muA']/phys['m']**2)
        if debug:
            print(Grid.shape)
            print(J.shape)

        deltaB = bs.deltaB('deltaB', x0, Grid, J, V_char = dx*dy*dz)
        return deltaB
        #return bs.B_EW(X0, Grid, J, Npole, dx*dy*dz)

    if print_output:
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

    if print_output:
        tf = t_module.time()
    #print(B_slices.shape)

    if print_output:
        print('Nx, Ny, Nz = {0:d}, {1:d}, {2:d}'.format(Nx,Ny,Nz))
        print('fullVolume = ' + str(fullVolume))
        print('fineVolume = ' + str(fineVolume))
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

        string = '\n X[0], X[-1], dx; Y[0], Y[-1], dy; Z[0], Z[-1], dz; Btot[0], Btot[1], Btot[2], np.linalg.norm(Btot) = {0:f}, {1:f}, {2:f}, {3:f}, {4:f}, {5:f}, {6:f}, {7:f}, {8:f}, {9:f}, {10:f}, {11:f}, {12:f}'.format(X[0], X[-1], dx, Y[0], Y[-1], dy, Z[0], Z[-1], dz, Btot[0], Btot[1], Btot[2], np.linalg.norm(Btot))
        print(string)
        if tpad(time) == (2003, 11, 20, 7, 0, 0, 0) and (mlat, mlon) == (57.50, 176.00):
            datafname = conf['run_path_derived'] + 'biot_savart_kameleon_data_2003:11:20T07:00:00.txt'
            f = open(datafname,'a') # append only mode
            f.write(string)
            f.close()

    return Btot
