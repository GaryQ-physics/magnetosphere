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
Nx_main = 3*27
Nx_tail = 30

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

def ex_data(kam,interp, variable, x,y,z, X0,Npole):
    if np.sqrt(x**2+y**2+z**2)<1e-4: return 0.
    # Get data from file, interpolate to point
    if('dB' in variable):
        J=np.array([ex_data(kam,interp, 'jx',x,y,z, X0,Npole),ex_data(kam,interp, 'jy',x,y,z, X0,Npole),ex_data(kam,interp, 'jz',x,y,z, X0,Npole)])
        R=np.array([x,y,z])-X0
        B=np.cross(J,R)/(np.linalg.norm(R)**3)
        if(variable=='dB_dV'): return np.linalg.norm(B)
        a2=np.cross(Npole,X0)
        if(variable=='dBlon_dV'): return np.dot(B,a2)/np.linalg.norm(a2)
        a1=np.cross(X0,a2)
        if(variable=='dBlat_dV'): return np.dot(B,a1)/np.linalg.norm(a1)
    kam.loadVariable(variable)
    data = interp.interpolate(variable, x, y, z)
    return data

def Compute(Event, var):
    year,month,day,hours,minutes,seconds,MLONdeg,MLATdeg = Event
    Time = [year,month,day,hours,minutes,seconds]
    MLON = MLONdeg*deg
    MLAT = MLATdeg*deg
    X0 = ps.MAGtoGSM([1.,MLATdeg,MLONdeg],Time,'sph','car')
    Npole = ps.GEOtoGSM([0.,0.,1.],Time,'car','car')

    # datafile
    filename = conf["f_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % (year,month,day,hours,minutes,seconds) + '.out.cdf'

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

    B2, B3, B1 = np.meshgrid(Y,Z,X)
    #B1, B2, B3 = np.meshgrid(X,Y,Z)

    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')
    B3 = B3.flatten(order='C')
    B = np.column_stack((B1,B2,B3))
    A = (np.nan)*np.empty((B1.size,))
    for l in range(A.size):
        A[l] = ex_data(kameleon,interpolator, var, B[l,0], B[l,1], B[l,2], X0,Npole)

    # close kameleon ---------------------
    kameleon.close()
    print("Closed " + filename)
    #-------------------------------
    return [A,B]

def writevtk(Event, var):
    A,B = Compute(Event,var)
    fname = conf["m_path"] + 'magnetosphere/data/kameleon_structured_grid_' + var + '.vtk'

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
    f.write('SCALARS point_scalars float 1\n')
    f.write('LOOKUP_TABLE default\n')
    np.savetxt(f, A)
    f.close()
    print("Wrote " + fname)
