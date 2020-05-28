# kameleon_structured_grid_write

import sys
import os
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')
sys.path.append(conf["m_path"] + 'magnetosphere/events/')
import _CCMC as ccmc
import pos_sun as ps

# run parameters
Ny = 3*30
Nz = 3*30
Nx_main = 2*3*27
Nx_tail = 20

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
kg = Tesla*A*s**2

mu0 = 1.2566370614e-6*kg*m/((s**2)*(A**2))

def ex_data(kam,interp, variable, x,y,z, X0,Npole, V_char = 1.):
    if np.sqrt(x**2+y**2+z**2)<1e-4: return 0.
    # Get data from file, interpolate to point
    if('dB' in variable):
        if np.sqrt(x**2+y**2+z**2)<1.5: return 0.
        J = np.array([ex_data(kam, interp, 'jx', x, y, z, X0, Npole), ex_data(kam, interp, 'jy', x, y, z, X0, Npole), ex_data(kam, interp, 'jz', x, y, z, X0, Npole)])
        J = J*(muA/m**2)
        R = np.array([x, y, z])-X0
        R = R*m
        B = (mu0/(4*np.pi))*np.cross(J, R)/(np.linalg.norm(R)**3)
        dBnT = B*V_char/(nT)
        if(variable=='dB_dV'): return np.linalg.norm(dBnT)
        a2 = np.cross(Npole, X0)
        if(variable=='dBlon_dV'): return np.dot(dBnT, a2)/np.linalg.norm(a2)
        a1 = np.cross(X0, a2)
        if(variable=='dBlat_dV'): return np.dot(dBnT, a1)/np.linalg.norm(a1)
    kam.loadVariable(variable)
    data = interp.interpolate(variable, x, y, z)
    return data

def Compute(Event, var):
    year,month,day,hours,minutes,seconds,MLONdeg,MLATdeg = Event
    time = [year,month,day,hours,minutes,seconds]
    MLON = MLONdeg*deg
    MLAT = MLATdeg*deg
    X0 = ps.MAGtoGSM([1.,MLATdeg,MLONdeg],time,'sph','car')
    Npole = ps.GEOtoGSM([0.,0.,1.],time,'car','car')

    # datafile
    filename = conf["run_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % (year,month,day,hours,minutes,seconds) + '.out.cdf'

    # open kameleon ---------------
    kameleon = ccmc.Kameleon()
    kameleon.open(filename)
    print(filename, "Opened " + filename)
    interpolator = kameleon.createNewInterpolator()
    #-----------------------------

    X = np.concatenate((np.linspace(-200,-20.05,Nx_tail), np.linspace(-20.,15.,Nx_main) ))
    Nx = X.size
    print('Nx= ',Nx)
    Y = np.linspace(-10.,10.,Ny)
    Z = np.linspace(-10.,10.,Nz)
    dV_grid = (X[-1] - X[-2])*(Y[1] - Y[0])*(Z[1] - Z[0])*R_e**3

    B2, B3, B1 = np.meshgrid(Y,Z,X)
    #B1, B2, B3 = np.meshgrid(X,Y,Z)

    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')
    B3 = B3.flatten(order='C')
    B = np.column_stack((B1,B2,B3))
    Aa = (np.nan)*np.empty((B1.size,))
    for l in range(Aa.size):
        Aa[l] = ex_data(kameleon,interpolator, var, B[l,0], B[l,1], B[l,2], X0,Npole, V_char = dV_grid)

    # close kameleon ---------------------
    kameleon.close()
    print("Closed " + filename)
    #-------------------------------
    return [Aa,B]

def writevtk(Event, var):
    A,B = Compute(Event,var)
    year,month,day,hours,minutes,seconds,MLONdeg,MLATdeg = Event
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d' % (year,month,day,hours,minutes,seconds)
    fname = conf["run_path_derived"] + 'kameleon_structured_grid_' + var + tag + '.vtk'

    Nx=Nx_main+Nx_tail
    print('also Nx= ',Nx)

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
    np.savetxt(f, A)
    f.close()
    print("Wrote " + fname)
