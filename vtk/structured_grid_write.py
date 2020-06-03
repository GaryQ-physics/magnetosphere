import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import _CCMC as ccmc
import pos_sun as ps

# run parameters
'''
Ny = 3*30
Nz = 3*30
Nx_main = 2*3*27
Nx_tail = 20
'''
dx = 0.2
dy = 0.2
dz = 0.2
dx_tail = 10.

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

def ex_data(kam,interp, variable, x,y,z, X0, Npole, V_char = 1.):
    if np.sqrt(x**2+y**2+z**2)<1e-4: return 0.
    # Get data from file, interpolate to point
    if('dB' in variable):
        if np.sqrt(x**2+y**2+z**2)<1.5: return 0.
        J = np.array([ex_data(kam, interp, 'jx', x, y, z, X0, Npole), 
                      ex_data(kam, interp, 'jy', x, y, z, X0, Npole), 
                      ex_data(kam, interp, 'jz', x, y, z, X0, Npole)])
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

def Compute(Event, var, calcTotal=False):
    total = 0.
    #year,month,day,hours,minutes,seconds,MLONdeg,MLATdeg = Event
    time = Event[0:6]
    MLON = Event[6]
    MLAT = Event[7]

    X0 = ps.MAGtoGSM([1.,MLAT,MLON],time,'sph','car')
    Npole = ps.GEOtoGSM([0.,0.,1.],time,'car','car')

    # datafile
    filename = conf["run_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % tuple(time) + '.out.cdf'

    # open kameleon ---------------
    kameleon = ccmc.Kameleon()
    kameleon.open(filename)
    print(filename, "Opened " + filename)
    interpolator = kameleon.createNewInterpolator()
    #-----------------------------

    X = np.concatenate((np.arange(-200,-20.05,dx_tail), np.arange(-20.,15.,dx) ))
    Nx = X.size
    print('Nx = ',Nx)
    Y = np.arange(-10., 10., dy)
    Ny = Y.size
    Z = np.arange(-10., 10., dz)
    Nz = Z.size

    B2, B3, B1 = np.meshgrid(Y, Z, X)
    #B1, B2, B3 = np.meshgrid(X, Y, Z)

    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')
    B3 = B3.flatten(order='C')
    B = np.column_stack((B1, B2, B3))
    Aa = (np.nan)*np.empty((B1.size, ))
    for l in range(Aa.size):
        Aa[l] = ex_data(kameleon, interpolator, var, B[l, 0], B[l, 1], B[l, 2], X0, Npole, V_char = dx*dy*dz) # dx*dy*dz*R_e**3
        if calcTotal:
            total = total + Aa[l]
    # close kameleon ---------------------
    kameleon.close()
    print("Closed " + filename)
    #-------------------------------
    if calcTotal:
        print('total = ', total)
    return [Aa, B, Nx, Ny, Nz]

def writevtk(Event, var, calcTotal=False):
    Aa, B, Nx, Ny, Nz = Compute(Event, var, calcTotal=calcTotal)
    time = Event[0:6]
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d' % tuple(time)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    fname = conf["run_path_derived"] + subdir + 'structured_grid_' + var + tag + '.vtk'
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)

    print('also Nx = ',Nx)

    f = open(fname,'w')
    print("Writing " + fname)
    f.write('# vtk DataFile Version 3.0\n')
    f.write('Structured Grid ' + var + '\n')
    f.write('ASCII\n')
    f.write('DATASET STRUCTURED_GRID\n')
    f.write('DIMENSIONS ' + str(Nx) + ' ' + str(Ny) + ' ' + str(Nz) + '\n' )
    f.write('POINTS '+str(Nx*Ny*Nz)+' float\n')
    np.savetxt(f, B)
    f.write('\n')
    f.write('POINT_DATA ' + str(Nx*Ny*Nz) + '\n')
    f.write('SCALARS '+ var +' float 1\n')
    f.write('LOOKUP_TABLE default\n')
    np.savetxt(f, Aa)
    f.close()
    print("Wrote " + fname)
