# field_line_vtk_script

import sys
import os
import numpy as np

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
import pos_sun as ps
import cut_plane

from scipy.integrate import odeint
import _CCMC as ccmc

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
        IC.append(ps.MAGtoGSM([R, MLAT-delta, MLON], time, 'sph', 'car'))

    if debug:
        print(IC)

    # Trace field lines
    s_grid = np.linspace(0., 200., 2000)

    if debug:
        print(s_grid)
        print(s_grid.size)

    solns = (np.nan)*np.empty((s_grid.size, 3, len(IC)))
    for i in range(len(IC)):
        sol = odeint(cut_plane.dXds, IC[i], s_grid, args = (kameleon, interpolator))
        solns[:, :, i] = sol

    if debug:
        print('IC', IC)

    # restrict the field lines to stop when reaching 1*R_E from the origin
    solns_restr = [] # initialize list of np_arrays, one for each restricted field line
    for i in range(Nb+1):  # loop over field lines
        # define condition on the field line points
        #tr = np.logical_and(solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 >=1.,solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 <= 224.**2+128.**2+128.**2)
        # create the arrays of the restricted field line componentwise
        went_out = 0
        end_val = solns.shape[0]-1
        if debug:
            print end_val
        for k in range(solns.shape[0]):
            if solns[k, 0, i]**2 + solns[k, 1, i]**2 + solns[k, 2, i]**2 > 1.1**2: went_out = 1
            if solns[k, 0, i]**2 + solns[k, 1, i]**2 + solns[k, 2, i]**2 < 1.1**2 and went_out==1:
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
        if debug and i==Nb+1: print(solns[:, :, i])
        if debug and i==Nb+1: print(sol)

    kameleon.close()
    if debug:
        print solns_restr
    return solns_restr


def writevtk(Event):
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    time = Event[0:6]
    T = tuple(time)

    tag = '_%04d:%02d:%02dT%02d:%02d:%02d' % T

    solns_restr = Compute(Event)
    mlong_array = [0., 10., -10., 20., -20.] # per deg
    for i in range(Nb+1+Nlong):
        out_fname = conf["run_path_derived"] + 'field_line' + str(i) + tag + '.vtk'
        if i > Nb:
            mlong = mlong_array[i-Nb-1]
            mlat = np.linspace(-90,90,100)
            a = np.column_stack([np.ones(100, ), mlat, mlong*np.ones(100, )])
            sol = ps.GEOtoGSM(a, time, 'sph', 'car')
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
                f.write('%e %e %e\n' % (sol[k, 0], sol[k, 1], sol[k, 2]))
            f.write('LINES ' + '1' + ' ' + str(sol.shape[0] + 1) + '\n' )
            f.write(str(sol.shape[0]) + '\n')
            for k in range(sol.shape[0]):
                f.write(str(k) + '\n')
            f.close()
            if debug:
                print('Wrote ' + out_fname)
        else:
            from_list = solns_restr[i]
            sol = np.array(from_list)
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
