import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

from scipy.integrate import odeint
import cxtransform as cx
from scipy.interpolate import RegularGridInterpolator

from util import maketag
from probe import probe


xlims = [-100., 15.]
ylims = [-10., 10.]
zlims = [-15., 15.]
dx = 0.3
dy = 0.3
dz = 0.3


def data2d(time, parameter, X, Y, U, debug=False):

    Z = np.zeros(X.shape)

    # grid of the corresponding values of variable. To be color plotted
    Z = data_in_U(time, parameter,
                        X.flatten(), Y.flatten(), U)

    return Z.reshape(X.shape)


def data_in_U(time, variable, u, v, U):
    """Data in U coordinates"""
    
    def norm(x):
        return x/np.sqrt(np.sum(x))

    # Normalize the U vectors.
    U1 = U[0]
    U2 = U[1]
    U3 = U[2]
    if isinstance(u, float):
        X = u*U1 + v*U2
    else:
        X = np.einsum('i,j->ij', u, U1) + np.einsum('i,j->ij', v, U2)
    

    if 'b' in variable:
        B = probe(time, X, var=['bx', 'by', 'bz'])
        if variable == 'bu1':
            return np.dot(B, U1)
        if variable == 'bu2':
            return np.dot(B, U2)
        if variable == 'bu3':
            return np.dot(B, U3)
    else:
        return probe(time, X, var=variable)

def dXds(X, s, sign, Bx_interp, By_interp, Bz_interp):

    if xlims[0]<X[0]<xlims[1] and ylims[0]<X[1]<ylims[1] and zlims[0]<X[2]<zlims[1]:
        #print('IN')
        B = np.array([Bx_interp(X)[0], By_interp(X)[0], Bz_interp(X)[0]])
        Bm = np.linalg.norm(B)
        if 1e-9 < Bm < 1e+7:
            return (sign/Bm)*B
    return [0., 0., 0.]


def fieldlines(time, mag, fieldvar='b', s_grid=None, max_iterations=100, debug=False):
    """Trace field lines from start points

    Parameters
    ----------
    time : list 
        or other type that is accepted by probe (cannot be 2d array)
    mag : numpy.array or list
        (Nx3) spherical mag coordinates of N points, which are the start points
        for N field lines
    s_grid : TYPE, optional
        DESCRIPTION. The default is np.arange(0., 200., 0.1)
    debug : TYPE, optional, boolean
        DESCRIPTION. The default is False.
    max_iterations : TYPE, optional, int
        stops integration after max_iterations amount, even if
        line goes out of bounds. The default is 100

    Returns
    -------
    ret : list of N arrays
        N field lines 

    """

    mag = np.array(mag)
    if len(mag.shape) == 1:
        mag = [mag]

    # Trace 3-D field line
    # TODO: Consider using
    # https://github.com/spacepy/spacepy/blob/master/spacepy/pybats/trace2d.py
    
    
    # make grid and interpolator on grid using probe ###################
    no_origin = xlims[0] > 0. or xlims[1] < 0. or ylims[0] > 0. or ylims[1] < 0. or zlims[0] > 0. or zlims[1] < 0.
    if no_origin:
        print('WARNING: grid does not contain origin')
        X = np.arange(xlims[0], xlims[1]+dx, dx)
        Y = np.arange(ylims[0], ylims[1]+dy, dy)
        Z = np.arange(zlims[0], zlims[1]+dz, dz)
    else:
        X = np.concatenate([ -np.flip(np.delete(np.arange(0., -xlims[0]+dx, dx), 0), 0) , np.arange(0., xlims[1]+dx, dx) ])
        Y = np.concatenate([ -np.flip(np.delete(np.arange(0., -ylims[0]+dy, dy), 0), 0) , np.arange(0., ylims[1]+dy, dy) ])
        Z = np.concatenate([ -np.flip(np.delete(np.arange(0., -zlims[0]+dz, dz), 0), 0) , np.arange(0., zlims[1]+dz, dz) ])
    Nx = X.size
    Ny = Y.size
    Nz = Z.size
    
    
    G2, G1, G3 = np.meshgrid(Y, X, Z) # different than in make_grid
    P = np.column_stack( (G1.flatten(), G2.flatten(), G3.flatten()) )
    Fx = probe(time, P, var=fieldvar+'x').reshape(Nx, Ny, Nz)
    Fy = probe(time, P, var=fieldvar+'y').reshape(Nx, Ny, Nz)
    Fz = probe(time, P, var=fieldvar+'z').reshape(Nx, Ny, Nz)

    # https://stackoverflow.com/questions/21836067/interpolate-3d-volume-with-numpy-and-or-scipy
    Fx_interp = RegularGridInterpolator((X,Y,Z), Fx)
    Fy_interp = RegularGridInterpolator((X,Y,Z), Fy)
    Fz_interp = RegularGridInterpolator((X,Y,Z), Fz)
    ####################################################################
  
    
    if s_grid is None:
        s_grid = np.arange(0., 10., 0.1)
    # Trace field line for a total length of smax, and check if stop conditions
    # satified. If not satified, trace for another total length of smax.
    # Note that Python 3 version of integration library has stop function
    # that can be passed so this won't be needed.

    IC = cx.MAGtoGSM(mag, time[0:6], 'sph', 'car')
    ret = []
    linenum = 0
    for X0 in IC:
        if debug:
            print('linenum = ' + str(linenum))
        done = False
        solns = np.empty((0, 3)) # Combined solutions
        i = 0
        while not done:
            if debug:
                print('i = ' + str(i))
            soln = odeint(dXds, X0, s_grid, args=(-1, Fx_interp, Fy_interp, Fz_interp))
            R = soln[:, 0]**2+soln[:, 1]**2 + soln[:, 2]**2
            # define condition on the field line points
            # Find first location where soln steps out-of-bounds
            #tr = np.where( False == (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.) )        
            # Boolean array.


            tr = (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.)
            # RuntimeWarning: invalid value encountered in greater_equal


            # Indices where stop conditions satisfied
            tr_out = np.where(tr == False)
            if debug:
                print(tr)
            if tr_out[0].size > 0:
                # Stop condition found at least once. Use solution up to that point.s
                solns = np.vstack((solns, soln[0:tr_out[0][0] + 1, :]))
                done = True
            elif max_iterations == i + 1:
                solns = np.vstack((solns, soln))   # return soln   faster?
                done = True
            else:
                # New initial condition is stop point.
                X0 = soln[-1, :]
                # Append solution but exclude last value, which is the
                # new initial condition.
                solns = np.vstack((solns, soln[0:-1, :]))
            i = i + 1
        ret.append(solns)
        linenum += 1
            
    return ret


def unitvector(time, mag, debug=False):

    # Points where field line is returned on are points a distance
    # 0, 0.5, 1.0 along the field line from the starting point.
    s_grid = np.array([0., 0.5, 1.])
    # Compute (x, y, z) of points at s_grid values.
    sols = fieldlines(time, mag, s_grid=s_grid, debug=debug, max_iterations=1)
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
                'cut_plane_info_%.2f_%.2f' % (mlat, mlon) + maketag(time) + '.txt'
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
