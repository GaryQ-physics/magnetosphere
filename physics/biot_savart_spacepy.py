import sys
import os
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import numpy as np

from units_and_constants import phys
import cxtransform as cx
import biot_savart as bs

import spacepy.pybats.bats as bats


def run(data3d, time, mlat, mlon, debug=False, returnX=False, quick=False):
    """

    Returns
    -------
    (3,) numpy array
        Using data from a SWMF IDL-Binary file, returns the result of biot savart 
        integration on fine near earth grid given said SWMF file, for the field
        at x0 corresponding to mlat, mlon on earth.

    Parameters
    ----------
    data3d : spacepy.pybats.bats.Bats2d
        this is the spacepy object with the loaded SWMF IDL-Binary file where 
        the current values will be pulled from.
    time : tuple/list/array
        the date-time in y, m, d, h, m, (s, milisec) implicitly assumed to be 
        consistent with the time in the file loaded in data3d
    mlat : float
        magentic latitude of x0 (radius is 1 R_e).
    mlon : float
        magentic latitude of x0 
    debug : boolean, optional
        DESCRIPTION. The default is False.
    returnX : boolean, optional
        DESCRIPTION. 
        if true, insted of returning the result of integration, returns (N,3) 
        array of all the gridpoints in the fine grid.
        The default is False.
    quick : boolean, optional
        if True, instead of integrating on fine mesh, it quickens the computation
        by selecting all points inside 3d box of the same size. This pickes up
        a few extra points.
        The default is False.

    """
    # get the cell coordinates
    x = data3d['x']
    y = data3d['y']
    z = data3d['z']
    jx = data3d['jx']
    jy = data3d['jy']
    jz = data3d['jz']

    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    z = np.array(z, dtype=float)
    jx = np.array(jx, dtype=float)
    jy = np.array(jy, dtype=float)
    jz = np.array(jz, dtype=float)

    D = 3.96875
    #xax = np.linspace(-D, D, 128)
    #yax = np.linspace(-D, D, 128)
    #zax = np.linspace(-D, D, 128)
    dx = 0.0625
    dV = dx**3
    assert((D+D)/(128-1) == dx)

    if quick:
        Tr = np.all([-D<=x, x<=D, -D<=y, y<=D, -D<=z, z<=D], axis=0)
        J_sp = np.column_stack([jx[Tr], jy[Tr], jz[Tr]])
        X = np.column_stack([x[Tr], y[Tr], z[Tr]])
    else:
        X = np.empty((0, 3)) # Combined slices in x
        J_sp = np.empty((0, 3)) # Combined slices in x
        for i in range(128):
            tr = x == i*0.0625 - D

            S = np.column_stack([x[tr], y[tr], z[tr]]) # each slice
            X = np.vstack((X, S))

            S = np.column_stack([jx[tr], jy[tr], jz[tr]])
            J_sp = np.vstack((J_sp, S))


    x0 = cx.MAGtoGSM([1., mlat, mlon], time, 'sph', 'car')

    out_spacepy = bs.deltaB('deltaB', x0, X, J_sp*(phys['muA']/phys['m']**2), V_char = dV)
    B_spacepy = np.linalg.norm(out_spacepy)

    if debug:
        print('---')
        #print(B_spacepy)
        print('input: {0:d},{1:d},{2:d},{3:d},{4:d};{5:.2f},{6:.2f}'.format\
                    (time[0],time[1],time[2],time[3],time[4],mlat,mlon))
        print('dBmhd = ' + str(B_spacepy))

    if returnX:
        return X
    return out_spacepy
