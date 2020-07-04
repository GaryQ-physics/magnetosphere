# J_field_lines_write

import sys
import os
import numpy as np

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import pos_sun as ps
#import cut_plane
from cut_plane import ex_data, dXds

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


def Compute(Event):
    #Event = [year, month, day, hours, minutes, seconds, miliseconds, MLONdeg, MLATdeg]
    time = Event[0:7]
    MLON = Event[7]
    MLAT = Event[8]
    T = tuple(time)

    filename = conf["run_path"] + '3d__var_3_e' \
                + '%04d%02d%02d-%02d%02d%02d-%03d' % T + '.out.cdf'

    kameleon = ccmc.Kameleon()
    if debug:
        print("Opening " + filename)
    kameleon.open(filename)
    if debug:
        print("Opened " + filename)
    interpolator = kameleon.createNewInterpolator()

    ang = np.linspace(-30.*deg, 30.*deg, 5)
    #dist = np.linspace(6.5, 10., 3)
    D = np.linspace(1., 10., 90)
    IC = []
    IC.append([7.1, 0., 0.])
    
    for j in range(ang.size):
        rose = False
        fell = False
        start = 7.
        end = 10.
        for i in range(D.size):
            P = ex_data(kameleon, interpolator, 'p', D[i]*np.cos(ang[j]), 0., D[i]*np.sin(ang[j]))
            Jy = ex_data(kameleon, interpolator, 'jy', D[i]*np.cos(ang[j]), 0., D[i]*np.sin(ang[j]))
            if P > 5. and False:
                rose = True
                start = D[i]
            if rose and P < 1. and False:
                fell = True
                end = D[i]
                break
        dist = np.linspace(start-0.3, end+0.3, 6)
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
    kameleon.loadVariable('jx')
    kameleon.loadVariable('jy')
    kameleon.loadVariable('jz')
    for i in range(len(IC)):
        sol = odeint(dXds, IC[i], s_grid, args = (kameleon, interpolator, 'j', -1))
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
    time = Event[0:7]
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d.%03d' % tuple(time)

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



'''
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
'''

