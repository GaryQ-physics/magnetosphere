import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
sys.path.append(os.path.dirname('/home/gary/magnetosphere/misc/' ))

import _CCMC as ccmc
import pos_sun as ps
from cut_plane import ex_data
import biot_savart as bs
# run parameters

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
m = R_e/6.3781e+6  # R_e == 6.3781e+6*m  (note m < 10^-6 Re)
#kg = Tesla*A*s**2


def J_kameleon(kam, interp, X):
    Jkunits = (np.nan)*np.empty(X.shape)
    for k in range(X.shape[0]):
        Jkunits[k,0] = ex_data(kam, interp, 'jx', X[k,0], X[k,1], X[k,2])
        Jkunits[k,1] = ex_data(kam, interp, 'jy', X[k,0], X[k,1], X[k,2])
        Jkunits[k,2] = ex_data(kam, interp, 'jz', X[k,0], X[k,1], X[k,2])
    return Jkunits

def Compute(Event, var, calcTotal=False):
    #Event = [year, month, day, hours, minutes, seconds, miliseconds, MLONdeg, MLATdeg]
    time = Event[0:7]
    MLON = Event[7]
    MLAT = Event[8]

    X0 = ps.MAGtoGSM([1.,MLAT,MLON],time[0:6],'sph','car')
    Npole = ps.GEOtoGSM([0.,0.,1.],time[0:6],'car','car')

    # datafile
    filename = conf["run_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-%03d' % tuple(time) + '.out.cdf'

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
    Bgrid = np.column_stack((B1, B2, B3))
    if 'dB' in var:
        unit_v = np.array([1., 0., 0.])
        if '_EW' in var:
            a2 = np.cross(Npole, X0)
            a1 = np.cross(X0, a2)
            a1 = a1/np.linalg.norm(a1)
            a2 = a2/np.linalg.norm(a2)
            unit_v = a2
        unit_v = np.repeat([unit_v], Bgrid.shape[0], axis=0)

        deltaBnT = bs.deltaB('dB', Bgrid, X0, J=J_kameleon(kameleon, interpolator, Bgrid)*(muA/m**2)) #/nT
        Aa = np.einsum('ij,ij->i', deltaBnT, unit_v) #https://stackoverflow.com/questions/15616742/vectorized-way-of-calculating-row-wise-dot-product-two-matrices-with-scipy

    else:
        Aa = (np.nan)*np.empty((B1.size, ))
        for l in range(Aa.size):
            Aa[l] = ex_data(kameleon, interpolator, var, Bgrid[l, 0], Bgrid[l, 1], Bgrid[l, 2], X0, Npole, V_char = dx*dy*dz) # dx*dy*dz*R_e**3

    # close kameleon ---------------------
    kameleon.close()
    print("Closed " + filename)
    #-------------------------------

    if calcTotal:
        print('total = ', np.sum(Aa))
    return [Aa, Bgrid, Nx, Ny, Nz]

def writevtk(Event, var, calcTotal=False):
    Aa, B, Nx, Ny, Nz = Compute(Event, var, calcTotal=calcTotal)
    time = Event[0:7]
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d.%03d' % tuple(time)
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
