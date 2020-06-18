import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import _CCMC as ccmc
from units_and_constants import phys
import pos_sun as ps
from cut_plane import ex_data
import biot_savart as bs
from biot_savart_demo1 import J_kunits

# run parameters

dx = .5
dy = .5 # 0.2
dz = .5
dx_tail = 10.  # dx_tail not used if using bs.make_grid
Test = False
debug = False
new = False

def J_kameleon(kam, interp, X):
    Jkunits = (np.nan)*np.empty(X.shape)
    for k in range(X.shape[0]):
        if np.dot(X[k,:], X[k,:]) >= 1.5**2:
            Jkunits[k,0] = ex_data(kam, interp, 'jx', X[k,0], X[k,1], X[k,2])
            Jkunits[k,1] = ex_data(kam, interp, 'jy', X[k,0], X[k,1], X[k,2])
            Jkunits[k,2] = ex_data(kam, interp, 'jz', X[k,0], X[k,1], X[k,2])
        else:
            Jkunits[k,0] = 0.
            Jkunits[k,1] = 0.
            Jkunits[k,2] = 0.
    return Jkunits

def Compute(Event, var, calcTotal=False):
    """
    Given an event time (and location), outputs an array for a grid of GSM coordinates
        covering the magnetosphere, and a corresponding array of the physical quantity
        var at those grid points.

    Note: the location part of Event is only used if var is a quantity measured
        relative to the event location, for example dB.
    """

    #Event = [year, month, day, hours, minutes, seconds, miliseconds, MLONdeg, MLATdeg]
    time = Event[0:7]
    MLON = Event[7]
    MLAT = Event[8]

    if Test:
        X0 = np.array([0., 0.75, 0.])
        Npole = np.array([0., 0., 1.])
    else:
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

    '''
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
    '''
    #tail = -200.
    tail = -100.

    ret = bs.make_grid([tail, 15.], [-10., 10.], [-10., 10.], dx, dy, dz)
    Xgrid = ret[0]
    
    print('X0=',X0)
    if 'dB' in var or var=='J':
        unit_v = np.array([1., 0., 0.])
        if '_EW' in var:
            a2 = np.cross(Npole, X0)
            a1 = np.cross(X0, a2)
            a1 = a1/np.linalg.norm(a1)
            a2 = a2/np.linalg.norm(a2)
            unit_v = a2
        unit_v = np.repeat([unit_v], Xgrid.shape[0], axis=0)
        
        if Test:
            Jin = J_kunits(Xgrid)*(phys['muA']/phys['m']**2)
        else:
            Jin = J_kameleon(kameleon, interpolator, Xgrid)*(phys['muA']/phys['m']**2)

        if debug:
            print('unit_v=',unit_v)
            print('Jin=',Jin)
            print('Jin sum=', np.sum(Jin))

        if var=='J':
            Aa = Jin
        else:
            deltaBnT = bs.deltaB('dB', X0, Xgrid, Jin, V_char = dx*dy*dz) #/nT
            Aa = np.einsum('ij,ij->i', deltaBnT, unit_v) #https://stackoverflow.com/questions/15616742/vectorized-way-of-calculating-row-wise-dot-product-two-matrices-with-scipy
        if debug:
            print('Aa=',Aa)
    else:
        Aa = (np.nan)*np.empty((Xgrid.shape[0], ))
        for l in range(Aa.size):
            Aa[l] = ex_data(kameleon, interpolator, var, Xgrid[l, 0], Xgrid[l, 1], Xgrid[l, 2])

    # close kameleon ---------------------
    kameleon.close()
    print("Closed " + filename)
    #-------------------------------

    if calcTotal:
        print('Test=',Test)
        print('dx,dy,dz = ', dx, dy, dz)
        print('total = ', np.sum(Aa))
        print(Xgrid.shape)
        print(Aa.shape)
    return [Aa, Xgrid, ret[1], ret[2], ret[3]]

def writevtk(Event, var, calcTotal=False):
    Aa, Bb, Nx, Ny, Nz = Compute(Event, var, calcTotal=calcTotal)
    time = Event[0:7]
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d.%03d' % tuple(time)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    fname = conf["run_path_derived"] + subdir + 'structured_grid_' + var + tag + '.vtk'
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)

    if debug:
        print('also Nx = ',Nx)
    if new:
        Xind = np.array([i for i in range(a)])
        Yind = np.array([j for j in range(b)])
        Zind = np.array([k for k in range(c)])
        C2, C3, C1 = np.meshgrid(Yind, Zind, Xind)
        C1 = C1.flatten(order='C')
        C2 = C2.flatten(order='C')
        C3 = C3.flatten(order='C')
        ind_p_np = np.column_stack((C1, C2, C3))

        X = np.array([1.5*i for i in range(a)])
        Y = np.array([1.5*j for j in range(b)])
        Z = np.array([1.5*k for k in range(c)])
        B2, B3, B1 = np.meshgrid(Y, Z, X)
        #B1, B2, B3 = np.meshgrid(X, Y, Z) # seems more natural but doesnt work with vtk structured_grid format
        B1 = B1.flatten(order='C')
        B2 = B2.flatten(order='C')
        B3 = B3.flatten(order='C')
        p_np = np.column_stack((B1, B2, B3))

        fvals = 0.1*p_np[:,0]*p_np[:,1]*p_np[:,2]

        Xc = np.array([i for i in range(a-1)])
        Yc = np.array([j for j in range(b-1)])
        Zc = np.array([k for k in range(c-1)])
        D2, D3, D1 = np.meshgrid(Yc, Zc, Xc)
        D1 = D1.flatten(order='C')
        D2 = D2.flatten(order='C')
        D3 = D3.flatten(order='C')
        cell_inds = np.column_stack((D1, D2, D3))

        shifts = np.array( [[0, 0, 0],
                            [1, 0, 0],
                            [1, 1, 0],
                            [0, 1, 0],
                            [0, 0, 1],
                            [1, 0, 1],
                            [1, 1, 1],
                            [0, 1, 1]] )

        cust_cells = np.nan*np.empty(((a-1)*(b-1)*(c-1), 8))
        for k in range((a-1)*(b-1)*(c-1)):
            for n in range(8):
                tup = tuple(cell_inds[k,:] + shifts[n,:])
                cust_cells[k,n] = np.where(np.all([ind_p_np[:,0]==tup[0],ind_p_np[:,1]==tup[1],ind_p_np[:,2]==tup[2]],axis=0))[0][0]

        c_np = [('hexahedron', cust_cells.astype(int))]
        fvals_dict= {'sample_scalars': fvals}
        cust_mesh = meshio.Mesh(p_np, c_np, point_data=fvals_dict)
        print("Writing " + fname + 'with meshio')
        meshio.vtk.write(fname, cust_mesh, binary=False)
        print("Wrote " + fname)

    else:
        f = open(fname,'w')
        print("Writing " + fname)
        f.write('# vtk DataFile Version 3.0\n')
        f.write('Structured Grid ' + var + '\n')
        f.write('ASCII\n')
        f.write('DATASET STRUCTURED_GRID\n')
        f.write('DIMENSIONS ' + str(Nx) + ' ' + str(Ny) + ' ' + str(Nz) + '\n' )
        f.write('POINTS '+str(Nx*Ny*Nz)+' float\n')
        #f.write(B.tobytes())
        np.savetxt(f, Bb)
        f.write('\n')
        f.write('POINT_DATA ' + str(Nx*Ny*Nz) + '\n')
        if var=='J':
            f.write('VECTORS ' + var + ' float\n')
        else:
            f.write('SCALARS ' + var + ' float 1\n')
            f.write('LOOKUP_TABLE default\n')
        np.savetxt(f, Aa)
        f.close()
        print("Wrote " + fname)
