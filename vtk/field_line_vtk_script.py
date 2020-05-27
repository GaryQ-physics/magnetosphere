# field_line_vtk_script

import sys
import os
import numpy as np
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')
sys.path.append(conf["m_path"] + 'magnetosphere/events/')
from scipy.integrate import odeint
import _CCMC as ccmc
import pos_sun as ps

# run parameters
Nlong = 5
Nb = 6
sign = -1  # changes sign of magnetic field used to trace the field lines
debug = False

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
    if x**2 + y**2 + z**2 >=1.:
        return data
    else:
        return 0.

def data_in_U(kam, interp, variable, u, v, U1, U2):
    # Get the data in the U coordinates (defined by the cut plane vectors U1 and U2)
    x,y,z = u*U1+v*U2
    B = np.array([ex_data(kam, interp, 'bx', x, y, z), ex_data(kam, interp, 'by', x, y, z), ex_data(kam, interp, 'bz', x, y, z)])
    if variable == 'bu1':
        return np.dot(B, U1)
    if variable == 'bu2':
        return np.dot(B, U2)
    if variable == 'bu3':
        return np.dot(B, U3)
    else:
        return ex_data(kam, interp, variable, x, y, z)

def dXds(X, s, kam, interp):
    '''
    derivative function for field line ODE
    dx/ds = Bx(x,y,z)/Bm
    dy/ds = By(x,y,z)/Bm
    dz/ds = Bz(x,y,z)/Bm
    
    X = [x, y, z]
    B = [Bx, By, Bz]
    Bm = sqrt(Bx**2 + By**2 + Bz**2)
    s=arclength    
    '''
    B = np.array([ex_data(kam, interp, 'bx', X[0], X[1], X[2]), ex_data(kam, interp, 'by', X[0], X[1], X[2]), ex_data(kam, interp, 'bz', X[0], X[1], X[2])])
    Bm = np.sqrt(np.dot(B, B))
    if 1e-9<Bm<1e+7:
        return (sign/Bm)*B
    else:
        if debug:
            if Bm >= 1e+7: print('FIELD TOO HIGH')
            if Bm <= 1e-7: print('FIELD TOO LOW')
        return [0., 0., 0.]

def Compute(Event):
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    time = Event[0:6]
    MLON = Event[6]
    MLAT = Event[7]
    T=tuple(time)

    filename = conf["f_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % T + '.out.cdf'

    # open kameleon ---------------
    kameleon = ccmc.Kameleon()
    kameleon.open(filename)
    print(filename, "Opened " + filename)
    interpolator = kameleon.createNewInterpolator()
    #-----------------------------

    R = 1.01
    eps = 3.
    IC = []
    for i in range(Nb+1):
        if i==0:
            delta = 0.
        elif i>Nb/2:
            #delta=(i-Nb/2)*eps
            delta = (-i)*eps
        else:
            #delta=(i-Nb/2)*eps
            delta = (-i)*eps
        IC.append(ps.MAGtoGSM([R,MLAT-delta,MLON],time,'sph','car'))
    if debug: print(IC)

    # Trace field lines
    s_grid = np.linspace(0., 200., 2000)
    if debug: print(s_grid)
    if debug: print(s_grid.size)
    solns = (np.nan)*np.empty((s_grid.size, 3, len(IC)))
    for i in range(len(IC)):
        sol = odeint(dXds, IC[i], s_grid, args = (kameleon, interpolator))
        solns[:, :, i] = sol
    if debug: print('IC', IC)

    # initialize vectors for defining field line cut plane
    v1 = (np.nan)*np.empty((3, ))
    v2 = (np.nan)*np.empty((3, ))
    v3 = (np.nan)*np.empty((3, ))
    U1 = (np.nan)*np.empty((3, ))
    U2 = (np.nan)*np.empty((3, ))
    U3 = (np.nan)*np.empty((3, ))

    # restrict the field lines to stop when reaching 1*R_E from the origin
    solns_restr = [] # initialize list of np_arrays, one for each restricted field line
    for i in range(Nb+1):  # loop over field lines
        # define condition on the field line points
        #tr = np.logical_and(solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 >=1.,solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 <= 224.**2+128.**2+128.**2)
        # create the arrays of the restricted field line componentwise
        went_out = 0
        came_back = 0
        end_val = solns.shape[0]-1
        for k in range(solns.shape[0]):
            if solns[k, 0, i]**2 + solns[k, 1, i]**2 + solns[k, 2, i]**2 > 1.1**2: went_out = 1
            if solns[k, 0, i]**2 + solns[k, 1, i]**2 + solns[k, 2, i]**2 and went_out==1:
                came_back = 1
                end_val = k
                break
        tr = np.arange(end_val+1)
        solx = solns[:, 0, i]
        solx = solx[tr]
        soly = solns[:, 1, i]
        soly = soly[tr]
        solz = solns[:, 2, i]
        solz = solz[tr]
        # reasemble and add to the list
        sol = np.column_stack([solx, soly, solz])
        solns_restr.append(sol)
        if debug and i==Nb+1: print(solns[:, :, i])
        if debug and i==Nb+1: print(sol)

    # close kameleon ---------------------
    kameleon.close()
    print("Closed " + filename)
    #-------------------------------
    return solns_restr


def writevtk(Event):
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    time = Event[0:6]
    MLON = Event[6]
    MLAT = Event[7]
    T=tuple(time)

    tag = '_%04d:%02d:%02dT%02d:%02d:%02d' % T

    solns_restr = Compute(Event)
    mlong_array = [0., 10., -10., 20., -20.] # per deg
    for i in range(Nb+1+Nlong):
        out_fname = conf["m_path"] + 'magnetosphere/data/' + 'field_line' + str(i) + tag + '.vtk'
        if i > Nb:
            mlong = mlong_array[i-Nb-1]
            mlat = np.linspace(-90,90,100)
            a = np.column_stack([np.ones(100, ), mlat, mlong*np.ones(100, )])
            sol = ps.GEOtoGSM(a, time, 'sph', 'car')
            print('writing ' + out_fname)
            f = open(out_fname, 'w')
            f.write('# vtk DataFile Version 3.0\n')
            f.write('A dataset with one polyline and no attributes\n')
            f.write('ASCII\n')
            f.write('\n')
            f.write('DATASET POLYDATA\n')
            f.write('POINTS ' + str(sol.shape[0]) + ' float\n')
            for k in range(sol.shape[0]):
                f.write('%e %e %e\n'%(sol[k, 0], sol[k, 1], sol[k, 2]))
            f.write('LINES ' + '1' + ' ' + str(sol.shape[0] + 1) + '\n' )
            f.write(str(sol.shape[0]) + '\n')
            for k in range(sol.shape[0]):
                f.write(str(k) + '\n')
            f.close()
            print('closed ' + out_fname)
        else:
            from_list = solns_restr[i]
            sol = np.array(from_list)
            print('writing ' + out_fname)
            f = open(out_fname, 'w')
            f.write('# vtk DataFile Version 3.0\n')
            f.write('A dataset with one polyline and no attributes\n')
            f.write('ASCII\n')
            f.write('\n')
            f.write('DATASET POLYDATA\n')
            f.write('POINTS ' + str(sol.shape[0]) + ' float\n')
            for k in range(sol.shape[0]):
                f.write('%e %e %e\n'%(sol[k, 0], sol[k, 1], sol[k, 2]))
            f.write('LINES ' + '1' + ' ' + str(sol.shape[0] + 1) + '\n' )
            f.write(str(sol.shape[0]) + '\n')
            for k in range(sol.shape[0]):
                f.write(str(k) + '\n')
            print('closed ' + out_fname)
            f.close()
