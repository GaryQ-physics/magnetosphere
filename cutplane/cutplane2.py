import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

from scipy.integrate import odeint
import cxtransform as cx
from scipy.interpolate import RegularGridInterpolator

import util
from probe import GetRunData
from fieldlines import fieldlines
from make_grid import make_grid, make_axes
from units_and_constants import phys
import biot_savart as bs
from biot_savart_kameleon_interpolated_grid import toMAGLocalComponents #TODO, get rid


def data_in_plane(run, time, variable, plane, **kwargs):
    """Data in U coordinates

    isinstance(plane, dict) is True
    ulims, vlims, e1, e2, e3  in plane.keys()

    plane['e3'] is a unit normal to the plane represented by (3,) array (in some background coordinates, say gsm
    plane['e1'] is unit vector in the positive "u" direction, with limits plane['ulims']
    plane['e2'] is unit vector in the positive "v" direction, with limits plane['vlims']
    """

    axes = make_axes(plane['ulims'], plane['vlims'], None, None,
                     dx=plane['du'], dy=plane['dv'])
    grid_in_plane = make_grid(axes)[:,:2]

    Npts = grid_in_plane.shape[0]
    u = grid_in_plane[:, 0]
    v = grid_in_plane[:, 1]
    U1 = plane['e1']
    U2 = plane['e2']
    U3 = plane['e3']

    if isinstance(u, float):
        #X = u*U1 + v*U2
        raise ValueError
    #grid_in_space = np.einsum('i,j->ij', u, U1) + np.einsum('i,j->ij', v, U2)
    grid_in_space = np.empty((Npts,3))
    for i in range(Npts): #TODO, array it
        grid_in_space[i,:] = grid_in_plane[i,0]*U1 + grid_in_plane[i,1]*U2

    '''
    for i in range(Npts):
        for j in range(3):
            X[i,j] = u[i]*U1[j] + v[i]*U2[j]
            X[i,j] = grid_in_plane[i,0]*U1[j] + grid_in_plane[i,1]*U2[j]

    for i in range(Npts):
        X[i,:] = grid_in_plane[i,0]*U1 + grid_in_plane[i,1]*U2

    X = np.einsum('ij,kj->ki', np.column_stack([U1,U2]), grid_in_plane)
    X = np.einsum('ji,kj->ki', np.array([U1,U2]), grid_in_plane)

    '''

    if 'bu' in variable:
        B = GetRunData(run,time, X, 'b')
        if variable == 'bu1':
            return np.dot(B, U1)
        if variable == 'bu2':
            return np.dot(B, U2)
        if variable == 'bu3':
            return np.dot(B, U3)
    elif 'dB' in variable:
        assert('mlat_dB' in kwargs.keys())
        assert('mlon_dB' in kwargs.keys())
        if 'V_char' in plane.keys():
            V_char = plane['V_char']
            if V_char is None:
                V_char = 1
        else:
            V_char = 1

        J = GetRunData(run, time, grid_in_space, 'j')*(phys['muA']/phys['m']**2)
        x0 = cx.MAGtoGSM(np.array([1.,kwargs['mlat_dB'],kwargs['mlon_dB']]), time, 'sph', 'car')
        dB = bs.deltaB('dB', x0, grid_in_space, J, V_char=V_char)

        if variable == 'dB_Magnitude':
            #print('dB.shape =' + str(dB.shape))
            return np.sqrt(np.einsum('ij,ij->i',dB, dB))
        else:
            dB_loc = toMAGLocalComponents(time, kwargs['mlat_dB'], kwargs['mlon_dB'], dB)
            if variable == 'dB_north':
                return dB_loc[:, 0]
            if variable == 'dB_east':
                return dB_loc[:, 1]
            if variable == 'dB_down':
                return dB_loc[:, 2]
    else:
# Since _CCMC used in file filemeta, in order to run plotting, cannot
# pass library='kameleonV' since then both kameleon and kameleonV
# library='kameleon' in the same shell.
        return GetRunData(run, time, grid_in_space, variable)


def unitvector(run, time, mag, debug=False): #TODO: revisit
    raise RuntimeWarning ('WARNING: need to revisit')
    # Points where field line is returned on are points a distance
    # 0, 0.5, 1.0 along the field line from the starting point.
    s_grid = np.array([0., 0.5, 1.])
    # Compute (x, y, z) of points at s_grid values.
    sols = fieldlines(run, time, mag, s_grid=s_grid, debug=debug, max_iterations=1)
    ret = []
    for sol in sols:
        if debug: print(sol.shape)
        # initialize vectors for defining field line cut plane
        v1 = (np.nan)*np.empty((3, ))
        v2 = (np.nan)*np.empty((3, ))
        v3 = (np.nan)*np.empty((3, ))
        U1 = (np.nan)*np.empty((3, ))
        U2 = (np.nan)*np.empty((3, ))
        U3 = (np.nan)*np.empty((3, ))

        # Three vectors in from origin to point on field line.
        v1 = sol[0, :]
        v2 = sol[2, :]
        v3 = sol[1, :]

        # Define cut plane coordinates based on field line 
        # (U3 is normal to the plane)
        U2 = (v1 - v2)/np.linalg.norm(v1-v2)
        U3 = np.cross(v3 - v1, U2)

        if np.linalg.norm(U3) < 1e-3:
            print("WARNING: close to straight line")
        U3 = U3/np.linalg.norm(U3)
        U1 = np.cross(U2, U3)
        ret.append([U1, U2, U3])
    return ret


def writedata(time, mlat, mlon, debug=False):
    """Write output of unitvector() to file
    
    Calling unitvector() from ParaView does not work, so write output to
    txt file.
    """

    # Compute centered dipole unit vector in GSM at given time
    Mdipole = cx.MAGtoGSM([0., 0., 1.], time[0:6], 'car', 'car')
    U1, U2, U3 = unitvector(time, np.array([1., mlat, mlon]))[0]

    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)

    out_fname = conf["run_path_derived"] + subdir + \
                'cut_plane_info_%.2f_%.2f' % (mlat, mlon) + util.maketag(time) + '.txt'
    f = open(out_fname, 'w')
    
    print('Writing ' + out_fname)
    f.write('%.7e %.7e %.7e\n' % (Mdipole[0], Mdipole[1], Mdipole[2]))
    f.write('%.7e %.7e %.7e\n' % (U1[0], U1[1], U1[2]))
    f.write('%.7e %.7e %.7e\n' % (U2[0], U2[1], U2[2]))
    f.write('%.7e %.7e %.7e\n' % (U3[0], U3[1], U3[2]))
    f.close()

    if debug:
        print('Wrote ' + out_fname)
        print(time)
        print('Mdipole = ', Mdipole)
        print('U1 = ', U1)
        print('U2 = ', U2)
        print('U3 = ', U3)
