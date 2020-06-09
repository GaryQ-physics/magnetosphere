# cut_plane

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

from scipy.integrate import odeint
import _CCMC as ccmc
import pos_sun as ps

# run parameters
sign = -1  # changes sign of magnetic field used to trace the field lines

# units
hr = 1.
muA = 1.
R_e = 1.
nT = 1.

deg = (np.pi/180.)
amin = deg/60.
minn = hr/60.
s = minn/60.

A = 1e+6*muA
Tesla = 1e+9*nT
m = R_e/6.3781e+6  # R_e == 6.3781e+6*m
#kg = Tesla*A*s**2

#mu0 = 1.2566370614e-6*kg*m/((s**2)*(A**2))
mu0 = 1.970237314e-10*nT*R_e/muA

def ex_data_full(kam, interp, variable, x, y, z, X0, Npole, V_char = 1.):
    if np.sqrt(x**2+y**2+z**2)<1e-4: return 0.
    # Get data from file, 'Interate to point
    if('dB' in variable):
        if np.sqrt(x**2+y**2+z**2)<1.5: return 0.
        J = np.array([ex_data_full(kam, interp, 'jx', x, y, z, X0, Npole), 
                      ex_data_full(kam, interp, 'jy', x, y, z, X0, Npole), 
                      ex_data_full(kam, interp, 'jz', x, y, z, X0, Npole)])
        J = J*(muA/m**2)
        R = X0 - np.array([x, y, z])
        #R = R*R_e
        #dB_dV = (mu0/(4*np.pi))*np.cross(J, R)/(np.linalg.norm(R)**3)
        #dBnT = dB_dV*V_char/(nT)
        dBnT = V_char*(mu0/(4*np.pi))*np.cross(J, R)/(np.linalg.norm(R)**3)
        if(variable=='dB'):
            return np.linalg.norm(dBnT)
        a2 = np.cross(Npole, X0)
        if(variable=='dB_EW'):
            return np.dot(dBnT, a2)/np.linalg.norm(a2) # east west direction (east positive)
        a1 = np.cross(X0, a2)
        if(variable=='dB_NS'):
            return np.dot(dBnT, a1)/np.linalg.norm(a1) # north south direction (north positive)
    kam.loadVariable(variable)
    data = interp.interpolate(variable, x, y, z)
    return data


def ex_data(kam, interp, variable, x, y, z):
    """Load data from file, interpolate to point"""

    if (x**2 + y**2 + z**2 >= 1.):
        return ex_data_full(kam,interp, variable, x, y, z, 0, 0)
    else:
        return 0.


def dXds(X, s, kam, interp, var):
    """Derivative function for field line ODE

    dx/ds = Bx(x,y,z)/Bm
    dy/ds = By(x,y,z)/Bm
    dz/ds = Bz(x,y,z)/Bm
    
    X = [x, y, z]
    B = [Bx, By, Bz]
    Bm = sqrt(Bx**2 + By**2 + Bz**2)
    s = arclength    
    """
    
    B = np.array([ex_data(kam, interp, var + 'x', X[0], X[1], X[2]), 
                  ex_data(kam, interp, var + 'y', X[0], X[1], X[2]), 
                  ex_data(kam, interp, var + 'z', X[0], X[1], X[2])])
    Bm = np.sqrt(np.dot(B, B))
    if 1e-9 < Bm < 1e+7:
        return (sign/Bm)*B
    else:
        return [0., 0., 0.] # TODO: Return np.nan?


def Compute(Event, ret_sol=False, r=1.01, debug=False):
    #Event = [year, month, day, hours, minutes, seconds, miliseconds, MLONdeg, MLATdeg]
    time = Event[0:7]
    MLON = Event[7]
    MLAT = Event[8]

    filename = conf["run_path"] + '3d__var_3_e' \
                + '%04d%02d%02d-%02d%02d%02d-%03d' % tuple(time) + '.out.cdf'

    X0 = ps.MAGtoGSM([r, MLAT, MLON], time[0:6], 'sph', 'car')

    if debug:
        print(filename, "Opening " + filename)
    kameleon = ccmc.Kameleon()
    kameleon.open(filename)
    if debug:
        print(filename, "Opened " + filename)
    interpolator = kameleon.createNewInterpolator()

    # Trace field line
    # TODO: Input to function should include ds and max length
    s_grid = np.linspace(0., 10., 101.)

    soln = odeint(dXds, X0, s_grid, args=(kameleon, interpolator, 'b'))
    if False:
        print('X0 = ')
        print(X0)
        print('X0 dot X0 = {0:.1e}'.format(np.dot(X0, X0)))
        print('soln = ')
        print(soln)

    # initialize vectors for defining field line cut plane
    v1 = (np.nan)*np.empty((3, ))
    v2 = (np.nan)*np.empty((3, ))
    v3 = (np.nan)*np.empty((3, ))
    U1 = (np.nan)*np.empty((3, ))
    U2 = (np.nan)*np.empty((3, ))
    U3 = (np.nan)*np.empty((3, ))
    Mdipole = (np.nan)*np.empty((3, ))

    # define condition on the field line points
    tr = np.logical_and(soln[:, 0]**2+soln[:, 1]**2 + soln[:, 2]**2 >=1.,
                        soln[:, 0]**2+soln[:, 1]**2 + soln[:,2]**2 < 20.)

    # create the arrays of the restricted field line componentwise
    sol = soln[tr, :]
    if debug:
        print('sol = ', sol)
        
    # define vects for plane of main field line
    if False:
        # Old method. Use start, middle, and end.
        # Won't work as desired when field line is not closed.
        v1 = sol[0,:]
        v2 = sol[-1,:]
        half = int(sol.shape[0]/2)
        v3 = sol[half,:]

    #import pdb; pdb.set_trace()
    # New method. Use first three field line points.
    v1 = sol[0, :]
    v2 = sol[2, :]
    v3 = sol[1, :]

    # define cut plane coordinates based on main field line 
    # (U3 is normal to the plane)
    U2 = (v1 - v2)/np.linalg.norm(v1-v2)
    U3 = np.cross(v3 - v1, U2)

    if np.linalg.norm(U3) < 1e-3:
        print("WARNING: close to straight line")
    U3 = U3/np.linalg.norm(U3)
    U1 = np.cross(U2, U3)

    kameleon.close()
    if debug:
        print("Closed " + filename)

    # Compute centered dipole vector in GSM at given time
    Mdipole = ps.MAGtoGSM([0., 0., 1.], time[0:6], 'car', 'car')

    if ret_sol: 
        return [Mdipole, U1, U2, U3, sol]
    else:
        return [Mdipole, U1, U2, U3]


def writedata(Event, debug=False):
    """Write output of compute() to file
    
    Calling compute() from ParaView does not work, so write output to file.
    """
    #year,month,day,hours,minutes,seconds,milisec,MLONdeg,MLATdeg = Event
    time = Event[0:7]
    Mdipole,U1,U2,U3 = Compute(Event)

    tag = '_%04d:%02d:%02dT%02d:%02d:%02d.%03d' % tuple(time)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)

    MLON=Event[6]
    MLAT=Event[7]

    out_fname = conf["run_path_derived"] + subdir + 'cut_plane_info_%.2f_%.2f' %(MLAT, MLON) + tag + '.txt'
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
