# cut_plane

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config_paths import config
conf = config()
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')

from scipy.integrate import odeint
import _CCMC as ccmc
import pos_sun as ps

# run parameters
debug = False
sign = -1  # changes sign of magnetic field used to trace the field lines

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

def ex_data(kam, interp, variable, x, y, z):
    # Get data from file, interpolate to point
    kam.loadVariable(variable)
    data = interp.interpolate(variable, x, y, z)
    if( x**2 + y**2 + z**2 >=1.):
        return data
    else:
        return 0.


def dXds(X, s, kam, interp):
    """
    derivative function for field line ODE
    dx/ds = Bx(x,y,z)/Bm
    dy/ds = By(x,y,z)/Bm
    dz/ds = Bz(x,y,z)/Bm
    
    X = [x, y, z]
    B = [Bx, By, Bz]
    Bm = sqrt(Bx**2 + By**2 + Bz**2)
    s=arclength    
    """
    
    B = np.array([ex_data(kam,interp, 'bx', X[0],X[1],X[2]), 
                  ex_data(kam,interp, 'by', X[0],X[1],X[2]), 
                  ex_data(kam,interp, 'bz', X[0],X[1],X[2])])
    Bm = np.sqrt(np.dot(B,B))
    if 1e-9 < Bm <1e+7:
        return (sign/Bm)*B
    else:
        if debug:
            if(Bm >= 1e+7): print('FIELD TOO HIGH')
            if(Bm <= 1e-7): print('FIELD TOO LOW')
        return [0., 0., 0.] # TODO: Return np.nan?


def Compute(Event):
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    time = Event[0:6]
    MLON = Event[6]
    MLAT = Event[7]
    T=tuple(time)

    filename = conf["run_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % T + '.out.cdf'

    R = 1.01
    X0 = ps.MAGtoGSM([R,MLAT,MLON],time,'sph','car')

    print(filename, "Opening " + filename)
    kameleon = ccmc.Kameleon()
    kameleon.open(filename)
    print(filename, "Opened " + filename)
    interpolator = kameleon.createNewInterpolator()

    # Trace field line
    s_grid = np.linspace(0., 10., 100.)
    soln = odeint(dXds, X0, s_grid, args=(kameleon, interpolator))
    if debug:
        print(X0)
        print(np.dot(X0,X0))
        print(soln)

    # initialize vectors for defining field line cut plane
    v1=(np.nan)*np.empty((3,))
    v2=(np.nan)*np.empty((3,))
    v3=(np.nan)*np.empty((3,))
    U1=(np.nan)*np.empty((3,))
    U2=(np.nan)*np.empty((3,))
    U3=(np.nan)*np.empty((3,))
    Mdipole=(np.nan)*np.empty((3,))

    # define condition on the field line points
    tr = np.logical_and(soln[:,0]**2+soln[:,1]**2+soln[:,2]**2 >=1., soln[:,0]**2+soln[:,1]**2+soln[:,2]**2 < 20.)

    # create the arrays of the restricted field line componentwise
    solx=soln[:,0]
    solx=solx[tr]
    soly=soln[:,1]
    soly=soly[tr]
    solz=soln[:,2]
    solz=solz[tr]

    # reassemble
    sol=np.column_stack([solx,soly,solz])
    print(sol)

    # define vects for plane of main field line
    v1 = sol[0,:]
    v2 = sol[-1,:]
    half = int(sol.shape[0]/2)
    v3 = sol[half,:]

    # define cut plane coordinates based on main field line 
    # (U3 is normal to the plane)
    U2 = (v1-v2)/np.linalg.norm(v1-v2)
    Mdipole = ps.MAGtoGSM([0.,0.,1.],time,'car','car')
    U3 = np.cross(v3-v1, U2)/np.linalg.norm(np.cross(v3-v1, U2))
    U1 = np.cross(U2, U3)

    kameleon.close()
    print("Closed " + filename)

    return [Mdipole,U1,U2,U3]


def writevtk(Event):
    """Write output of compute() to file
    
    Calliong compute() from ParaView does not work, so write output to file.
    """
    year,month,day,hours,minutes,seconds,MLONdeg,MLATdeg = Event
    time = [year,month,day,hours,minutes,seconds]
    Mdipole,U1,U2,U3 = Compute(Event)
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d' % (year,month,day,hours,minutes,seconds)

    out_fname = conf["run_path_derived"] + 'cut_plane_info' + tag + '.txt'
    f = open(out_fname, 'w')
    print('Writing ' + out_fname)
    f.write('%.7e %.7e %.7e\n' % (Mdipole[0], Mdipole[1], Mdipole[2]))
    f.write('%.7e %.7e %.7e\n' % (U1[0], U1[1], U1[2]))
    f.write('%.7e %.7e %.7e\n' % (U2[0], U2[1], U2[2]))
    f.write('%.7e %.7e %.7e\n' % (U3[0], U3[1], U3[2]))
    f.close()
    print('Wrote ' + out_fname)

    print(time)
    print('Mdipole = ', Mdipole)
    print('U1 = ', U1)
    print('U2 = ', U2)
    print('U3 = ', U3)
