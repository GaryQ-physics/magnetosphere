# J_field_lines_write

import sys
import os
import numpy as np

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import pos_sun as ps
#import cut_plane
from B_field_lines_write import ex_data

from scipy.integrate import odeint
import _CCMC as ccmc

# run parameters
debug = False

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.


def dXds(X, s, kam, interp):    
    J = np.array([ex_data(kam, interp, 'jx', X[0], X[1], X[2]), 
                  ex_data(kam, interp, 'jy', X[0], X[1], X[2]), 
                  ex_data(kam, interp, 'jz', X[0], X[1], X[2])])
    Jm = np.sqrt(np.dot(J, J))
    if 1e-9 < Jm <1e+7:
        return (1./Jm)*J
    else:
        return [0., 0., 0.] # TODO: Return np.nan?


def Compute(Event):
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    time = Event[0:6]
    MLON = Event[6]
    MLAT = Event[7]
    T = tuple(time)

    filename = conf["run_path"] + '3d__var_3_e' \
                + '%04d%02d%02d-%02d%02d%02d-000' % T + '.out.cdf'

    kameleon = ccmc.Kameleon()
    if debug:
        print("Opening " + filename)
    kameleon.open(filename)
    if debug:
        print("Opened " + filename)
    interpolator = kameleon.createNewInterpolator()

    ang = np.linspace(-1., 1., 5)
    #dist = np.linspace(6.5, 10., 3)
    D = np.linspace(1., 10., 90)
    IC = []
    for j in range(ang.size):
        rose = False
        fell = False
        start = 7.
        end = 10.
        for i in range(D.size):
            P = ex_data(kameleon, interpolator, 'p', D[i]*np.cos(ang[j]), 0., D[i]*np.sin(ang[j]))
            if P > 6.:
                rose = True
                start = D[i]
            if rose and P < 2.:
                fell = True
                end = D[i]
                break
        dist = np.linspace(start-0.3, end+0.3, 4)
        for k in range(dist.size):
            dist[k]*np.cos(ang[j]), 0., dist[k]*np.sin(ang[j])
            IC.append([dist[k]*np.cos(ang[j]), 0., dist[k]*np.sin(ang[j])])

    '''
    B1, B2 = np.meshgrid(dist,ang)
    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')

    if debug:
        print('B1=',B1)

    x = B1*np.cos(B2)
    z = B1*np.sin(B2)
    y = 0.*x

    XYZ = np.column_stack((x, y, z))

    IC = list(XYZ)
    

    if debug:
        print(XYZ)
    '''

    # Trace field lines
    s_grid = np.linspace(0., 10., 200)

    if debug:
        print(s_grid)
        print(s_grid.size)

    solns = (np.nan)*np.empty((s_grid.size, 3, len(IC)))
    for i in range(len(IC)):
        sol = odeint(dXds, IC[i], s_grid, args = (kameleon, interpolator))
        solns[:, :, i] = sol

    if True:
        print('IC', IC)

    # restrict the field lines to stop when reaching 1*R_E from the origin
    solns_restr = [] # initialize list of np_arrays, one for each restricted field line
    for i in range(len(IC)):  # loop over field lines
        # define condition on the field line points
        #tr = np.logical_and(solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 >=1.,solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 <= 224.**2+128.**2+128.**2)
        # create the arrays of the restricted field line componentwise
        end_val = solns.shape[0]-1
        if debug:
            print end_val
        for k in range(solns.shape[0]):
            if solns[k, 0, i]**2 + solns[k, 1, i]**2 + solns[k, 2, i]**2 < 1.1**2:
                end_val = k
                break
        if debug:
            print end_val
        tr = np.arange(end_val+1)
        solx = solns[:, 0, i]
        if debug:
            print solx
            print tr
        solx = solx[tr]

        if debug:
            print solx
        soly = solns[:, 1, i]
        soly = soly[tr]
        solz = solns[:, 2, i]
        solz = solz[tr]

        # reasemble and add to the list
        sol = np.column_stack([solx, soly, solz])
        solns_restr.append(sol)

    kameleon.close()
    if debug:
        print solns_restr
    return solns_restr


def writevtk(Event):
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    time = Event[0:6]
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d' % tuple(time)

    solns_restr = Compute(Event)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)
    for sol_list in solns_restr:
        sol = np.array(sol_list)

        out_fname = conf["run_path_derived"] + subdir + 'J_field_line_%.2f_%.2f_%.2f.vtk' %(sol[0, 0], sol[0, 1], sol[0, 2])

        if debug: 
            print('Writing ' + out_fname)

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
        if debug:
            print('Wrote ' + out_fname)
