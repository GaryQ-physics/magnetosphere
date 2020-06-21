# B_field_lines_write

import sys
import os
import numpy as np

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import pos_sun as ps
import biot_savart as bs
import cut_plane

from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import odeint
import _CCMC as ccmc

global_x_min = -224.
global_x_max = 32.
global_y_min = -128.
global_y_max = 128.
global_z_min = -128.
global_z_max = 128.

# run parameters
#Nlong = 5
#Nb = 6
sign = -1  # changes sign of magnetic field used to trace the field lines
debug = False

xlims = [-100., 15.]
ylims = [-10., 10.]
zlims = [-15., 15.]
dx = 0.3
dy = 0.3
dz = 0.3

Grid = bs.make_grid(xlims, ylims, zlims, dx, dy, dz)[0]
no_origin = xlims[0] > 0. or xlims[1] < 0. or ylims[0] > 0. or ylims[1] < 0. or zlims[0] > 0. or zlims[1] < 0.
if no_origin:
    print('WARNING: grid does not contain origin')
    X = np.arange(xlims[0], xlims[1]+dx, dx)
    Y = np.arange(ylims[0], ylims[1]+dy, dy)
    Z = np.arange(zlims[0], zlims[1]+dz, dz)
else:
    X = np.concatenate([ -np.flip(np.delete(np.arange(0., -xlims[0]+dx, dx), 0), 0) , np.arange(0., xlims[1]+dx, dx) ])
    Y = np.concatenate([ -np.flip(np.delete(np.arange(0., -ylims[0]+dy, dy), 0), 0) , np.arange(0., ylims[1]+dy, dy) ])
    Z = np.concatenate([ -np.flip(np.delete(np.arange(0., -zlims[0]+dz, dz), 0), 0) , np.arange(0., zlims[1]+dz, dz) ])
Nx = X.size
Ny = Y.size
Nz = Z.size

filename = conf['run_path'] + '3d__var_3_e20031120-070000-000.out.cdf'
kameleon = ccmc.Kameleon()
if debug:
    print("Opening " + filename)
kameleon.open(filename)
if debug:
    print("Opened " + filename)
interpolator = kameleon.createNewInterpolator()
# https://stackoverflow.com/questions/21836067/interpolate-3d-volume-with-numpy-and-or-scipy
Bx = np.nan*np.empty((Nx,Ny,Nz))
By = np.nan*np.empty((Nx,Ny,Nz))
Bz = np.nan*np.empty((Nx,Ny,Nz))
kameleon.loadVariable('bx')
kameleon.loadVariable('by')
kameleon.loadVariable('bz')
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            Bx[i,j,k] = interpolator.interpolate('bx', X[i], Y[j], Z[k])
            By[i,j,k] = interpolator.interpolate('by', X[i], Y[j], Z[k])
            Bz[i,j,k] = interpolator.interpolate('bz', X[i], Y[j], Z[k])
kameleon.close()
if debug:
    print("Closed " + filename)

Bx_interp = RegularGridInterpolator((X,Y,Z), Bx)
By_interp = RegularGridInterpolator((X,Y,Z), By)
Bz_interp = RegularGridInterpolator((X,Y,Z), Bz)

def dXds(X, s, sign):
    if xlims[0]<X[0]<xlims[1] and ylims[0]<X[1]<ylims[1] and zlims[0]<X[2]<zlims[1]:
        #print('IN')
        B = np.array([Bx_interp(X)[0], By_interp(X)[0], Bz_interp(X)[0]])
        Bm = np.linalg.norm(B)
        if 1e-9 < Bm < 1e+7:
            return (sign/Bm)*B
    return [0., 0., 0.] # TODO: Return np.nan?

#print(Bx_interp(np.array([1,2,3])))
#print(dXds(np.array([1,2,3]), 1))

def Compute(Event, Nb, use_grid=True):
    #Event = [year, month, day, hours, minutes, seconds, miliseconds, MLONdeg, MLATdeg]
    time = Event[0:7]
    MLON = Event[7]
    MLAT = Event[8]
    T = tuple(time)

    filename = conf["run_path"] + '3d__var_3_e' \
                + '%04d%02d%02d-%02d%02d%02d-%03d' % T + '.out.cdf'

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
        IC.append(ps.MAGtoGSM([R, MLAT-delta, MLON], time[0:6], 'sph', 'car'))

    if debug:
        print(IC)

    # Trace field lines
    s_grid = np.linspace(0., 200., 2000)

    if debug:
        print(s_grid)
        print(s_grid.size)

    solns = (np.nan)*np.empty((s_grid.size, 3, len(IC)))
    for i in range(len(IC)):
        if use_grid:
            sol = odeint(dXds, IC[i], s_grid, args = (-1,))
        else:
            kameleon = ccmc.Kameleon()
            if debug:
                print("Opening " + filename)
            kameleon.open(filename)
            if debug:
                print("Opened " + filename)
            interpolator = kameleon.createNewInterpolator()
            kameleon.loadVariable('bx')
            kameleon.loadVariable('by')
            kameleon.loadVariable('bz')
            sol = odeint(cut_plane.dXds, IC[i], s_grid, args = (kameleon, interpolator, 'b', -1))
            kameleon.close()

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
            print(end_val)
        for k in range(solns.shape[0]):
            if solns[k, 0, i]**2 + solns[k, 1, i]**2 + solns[k, 2, i]**2 > 1.1**2: went_out = 1
            if solns[k, 0, i]**2 + solns[k, 1, i]**2 + solns[k, 2, i]**2 < 1.1**2 and went_out==1:
                end_val = k
                break
        if debug:
            print(end_val)
        tr = np.arange(end_val+1)
        solx = solns[:, 0, i]
        if debug:
            print(solx)
            print(tr)
        solx = solx[tr]

        if debug:
            print(solx)
        soly = solns[:, 1, i]
        soly = soly[tr]
        solz = solns[:, 2, i]
        solz = solz[tr]

        # reasemble and add to the list
        sol = np.column_stack([solx, soly, solz])
        solns_restr.append(sol)
        if debug and i==Nb+1: print(solns[:, :, i])
        if debug and i==Nb+1: print(sol)

    if debug:
        print(solns_restr)
    return solns_restr


def writevtk(Event, Nb=6):
    #Event == [year, month, day, hours, minutes, seconds, milisec, MLONdeg, MLATdeg]
    time = Event[0:7]
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d.%03d' % tuple(time)

    solns_restr = Compute(Event, Nb)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)
    for i in range(Nb+1):
        from_list = solns_restr[i]
        sol = np.array(from_list)

        v = ps.GSMtoMAG([sol[0, 0], sol[0, 1], sol[0, 2]], time[0:6], 'car', 'sph')
        if i==0:
            out_fname = conf["run_path_derived"] + subdir + 'B_field_line_event_%.2f_%.2f.vtk' %(Event[7], Event[6]) #(MLAT, MLON)
        else:
            out_fname = conf["run_path_derived"] + subdir + 'B_field_line_%.2f_%.2f.vtk' %(v[1], v[2]) #(MLAT, MLON)

        if debug: 
            print('Writing ' + out_fname)
            if i==0:
                print('%.2f_%.2f' %(v[1], v[2]) == '%.2f_%.2f' %(Event[7], Event[6])) # should be True
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
