import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

from scipy.integrate import odeint
import _CCMC as ccmc
import cxtransform as cx

from util import time2filename


def data2d(time, parameter, X, Y, U, debug=False):

    filename = time2filename(time)

    kameleon = ccmc.Kameleon()
    if debug:
        print("Opening " + filename)
    kameleon.open(filename)
    if debug:
        print("Opened " + filename)
    interpolator = kameleon.createNewInterpolator()

    Z = np.zeros(X.shape)
    for i in range(X.shape[0]): # note this is y_1d.size (NOT x)
        for j in range(X.shape[1]): 
            # grid of the corresponding values of variable. To be color plotted
            Z[i, j] = data_in_U(kameleon, interpolator, parameter,
                                X[i, j], Y[i, j], U)

    kameleon.close()
    if debug:
        print("Closed " + filename + "\n")

    return Z


def data_in_U(kam, interp, variable, u, v, U):
    """Data in U coordinates"""
    
    def norm(x):
        return x/np.sqrt(np.sum(x))

    # Normalize the U vectors.
    U1 = U[0]
    U2 = U[1]
    U3 = U[2]
    
    x, y, z = u*U1 + v*U2
    B = np.array([ex_data(kam, interp, 'bx', x, y, z), 
                  ex_data(kam, interp, 'by', x, y, z), 
                  ex_data(kam, interp, 'bz', x, y, z)])
    if variable == 'bu1':
        return np.dot(B, U1)
    if variable == 'bu2':
        return np.dot(B, U2)
    if variable == 'bu3':
        return np.dot(B, U3)
    else:

        return ex_data(kam, interp, variable, x, y, z)


def ex_data(kam, interp, variable, x, y, z):
    """Load data from file, interpolate to point"""
    if (x**2 + y**2 + z**2 < 1.):
        return 0
    kam.loadVariable(variable)
    data = interp.interpolate(variable, x, y, z)
    return data


def dXds(X, s, kam, interp, var):
    """Derivative function for field line ODE

    dx/ds = Fx(x,y,z)/Fm
    dy/ds = Fy(x,y,z)/Fm
    dz/ds = Fz(x,y,z)/Fm
    
    X = [x, y, z]
    F = [Fx, Fy, Fz]
    Fm = sqrt(Fx**2 + Fy**2 + Fz**2)
    s = arclength

    F is magnetic field for          var = 'b'
    F is current density field for   var = 'j'
    """

    # run parameters
    sign = -1  # changes sign of magnetic field used to trace the field lines
        
    B = np.array([ex_data(kam, interp, var + 'x', X[0], X[1], X[2]), 
                  ex_data(kam, interp, var + 'y', X[0], X[1], X[2]), 
                  ex_data(kam, interp, var + 'z', X[0], X[1], X[2])])
    Bm = np.sqrt(np.dot(B, B))
    if 1e-9 < Bm < 1e+7:
        return (sign/Bm)*B
    else:
        return [0., 0., 0.] # TODO: Return np.nan?


def fieldlines(time, mag, s_grid=None, debug=False):

    # Trace 3-D field line
    # TODO: Consider using
    # https://github.com/spacepy/spacepy/blob/master/spacepy/pybats/trace2d.py

    filename = conf["run_path"] + '3d__var_3_e' \
                + '%04d%02d%02d-%02d%02d%02d-%03d' % tuple(time) + '.out.cdf'

    if debug:
        print(filename, "Opening " + filename)
    kameleon = ccmc.Kameleon()
    kameleon.open(filename)
    if debug:
        print(filename, "Opened " + filename)
    interpolator = kameleon.createNewInterpolator()

    if s_grid is None:
        s_grid = np.arange(0., 200., 0.1)
    # Trace field line for a total length of smax, and check if stop conditions
    # satified. If not satified, trace for another total length of smax.
    # Note that Python 3 version of integration library has stop function
    # that can be passed so this won't be needed.
    done = False
    solns = np.empty((0, 3)) # Combined solutions
    X0 = cx.MAGtoGSM(mag, time[0:6], 'sph', 'car')
    while not done:
        soln = odeint(dXds, X0, s_grid, args=(kameleon, interpolator, 'b'))
        R = soln[:, 0]**2+soln[:, 1]**2 + soln[:, 2]**2
        # define condition on the field line points
        # Find first location where soln steps out-of-bounds
        #tr = np.where( False == (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.) )        
        # Boolean array.
        tr = (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.)
        # Indices where stop conditions satisfied
        tr_out = np.where(tr == False)
        if tr_out[0].size > 0:
            # Stop condition found at least once. Use solution up to that point.s
            solns = np.vstack((solns, soln[0:tr_out[0][0] + 1, :]))
            done = True
        else:
            # New initial condition is stop point.
            X0 = soln[-1, :]
            # Append solution but exclude last value, which is the
            # new initial condition.
            solns = np.vstack((solns, soln[0:-1, :]))
    
    kameleon.close()
    
    return solns


def unitvector(time, mag, debug=False):

    # Points where field line is returned on are points a distance
    # 0, 0.5, 1.0 along the field line from the starting point.
    s_grid = np.array([0., 0.5, 1.])
    # Compute (x, y, z) of points at s_grid values.
    sol = fieldlines(time, mag, s_grid=s_grid, debug=debug)
    
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

    return [U1, U2, U3]


def writedata(time, mlat, mlon, debug=False):
    """Write output of unitvector() to file
    
    Calling unitvector() from ParaView does not work, so write output to
    txt file.
    """

    # Compute centered dipole unit vector in GSM at given time
    Mdipole = cx.MAGtoGSM([0., 0., 1.], time[0:6], 'car', 'car')

    Mdipole, U1, U2, U3 = unitvector(time)

    tag = '_%04d:%02d:%02dT%02d:%02d:%02d.%03d' % tuple(time)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)

    out_fname = conf["run_path_derived"] + subdir + \
                'cut_plane_info_%.2f_%.2f' % (mlat, mlon) + tag + '.txt'
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
