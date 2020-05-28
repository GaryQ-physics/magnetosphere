# cut_plane_python_plotting

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config_paths import config
conf = config()

import matplotlib.pyplot as plt
# Following needed for projection='3d', but PyFlakes will warn that not used
from scipy.integrate import odeint

sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')
import _CCMC as ccmc
import pos_sun as ps

from cut_plane import ex_data, dXds, Compute

def data_in_U(kam, interp, variable, u, v, U1, U2, U3):
    # Get the data in the U coordinates (defined by the cut plane vectors U1 and U2)
    x,y,z = u*U1 + v*U2
    B = np.array([ex_data(kam, interp, 'bx', x, y, z), 
                  ex_data(kam, interp, 'by', x, y, z), 
                  ex_data(kam, interp, 'bz', x, y, z)])
    if variable == 'bu1':
        return np.dot(B, U1)
    if variable == 'bu2':
        return np.dot(B, U2)
    if variable == 'bu3':
        return np.dot(B, U3)
    else:
        return ex_data(kam, interp, variable, x, y, z)


def plot(time, pos, plane_vs, parameter, xlim=[0,4], ylim=[-3,3], nx=50, ny=50, png=True):
    # plot(time, [r, mlat, mlong], None, 'p')
    # plot(time, [GSMx,GSMy,GSMz], [v1, v2], 'p')
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    
    # run parameters
    debug = False

    # Plot title
    title = 'SCARR5 ' + '%04d%02d%02dT%02d%02d%02d' % tuple(time)
    title = title + "\n" + "[mlat,mlon]=[{0:.1f}, {1:.1f}]".format(pos[1], pos[2])

    filename = '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % tuple(time)

    filename_in = conf["run_path"] + filename + '.out.cdf'
    filename_out = conf["run_path_derived"] + filename + '.png'
    
    r = 1.01
    X0 = ps.MAGtoGSM(pos, time, 'sph', 'car')

    # open kameleon
    kameleon = ccmc.Kameleon()
    if debug:
        print("Opening " + filename_in)
    kameleon.open(filename_in)
    if debug:
        print("Opened " + filename_in)
    interpolator = kameleon.createNewInterpolator()

    parameter_unit = kameleon.getVisUnit(parameter)
    parameter_unit = parameter_unit.replace('mu', '\\mu ')

    if plane_vs==None:
        '''
        Event = time + pos[1:3]
        vects = Compute(Event)
        U1 = vects[1]
        U2 = vects[2]
        U3 = vects[3]
        '''

        # Trace field line
        s_grid = np.linspace(0., 10., 100)
        soln = odeint(dXds, X0, s_grid, args=(kameleon, interpolator))
        if debug:
            print(X0)
            print(np.dot(X0, X0))
            print(soln)

        # initialize vectors for defining field line cut plane
        v1 = (np.nan)*np.empty((3, ))
        v2 = (np.nan)*np.empty((3, ))
        v3 = (np.nan)*np.empty((3, ))
        U1 = (np.nan)*np.empty((3, ))
        U2 = (np.nan)*np.empty((3, ))
        U3 = (np.nan)*np.empty((3, ))

        # define restriction condition on the field line points
        tr = np.logical_and(soln[:, 0]**2 + soln[:, 1]**2 + soln[:, 2]**2 >= 1.,
                            soln[:, 0]**2 + soln[:, 1]**2 + soln[:, 2]**2 < 20.)

        # restricted field lines
        sol = soln[tr, :]

        # define vects for plane of main field line
        v1 = sol[0, :]  # First point on field line
        v2 = sol[-1, :] # Last point on field line
        half = int(sol.shape[0]/2) 
        v3 = sol[half, :] # Approximate mid-point on field line

        # define cut plane coordinates based on main field line 
        # (U3 is normal to the plane)
        U2 = (v1-v2)/np.linalg.norm(v1-v2)
        U3 = np.cross(v3-v1, U2)/np.linalg.norm(np.cross(v3-v1, U2))
        U1 = np.cross(U2, U3)
    else:
        U1 = np.array(plane_vs[0])
        U2 = np.array(plane_vs[1])
        U3 = np.cross(U1, U2)

    x_1d = np.linspace(xlim[0], xlim[1], nx)
    y_1d = np.linspace(ylim[0], ylim[1], ny)
    X, Y = np.meshgrid(x_1d, y_1d) # grid of points on the cutplane
    Z = np.zeros((nx, ny))
    for i in range(nx):
        for j in range(ny):
            # grid of the corresponding values of variable. To be color plotted
            Z[i, j] = data_in_U(kameleon, interpolator, parameter, X[i, j], Y[i, j], U1, U2, U3)

    kameleon.close()
    if debug:
        print("Closed " + filename_in + "\n")

    # Plotting
    plt.clf()
    fig = plt.figure(dpi = 200)
    ax2 = fig.add_subplot(1, 2, 2)

    # Plot cut plane data
    ax2.set_title(title, fontsize=10)
    ax2.set(xlabel = "Tailward distance [$R_E$]")
    ax2.set(ylabel = "Northward distance [$R_E$]")
    ax2.axis('square')
    pcm = ax2.pcolormesh(X, Y, Z)

    # Reason for choice of fraction and pad:
    # https://stackoverflow.com/a/39948312/1491619
    cb = fig.colorbar(pcm, ax = ax2, fraction = 0.04, pad = 0.08)
    cb.ax.set_title(parameter + ' [$' + parameter_unit + '$]', ha='left')
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])

    # Add field line to 2D plot
    if plane_vs==None:
        solCut=np.zeros((sol.shape[0], 2))
        for k in range(sol.shape[0]):
            solCut[k, 0] = np.dot(sol[k, :], U1)
            solCut[k, 1] = np.dot(sol[k, :], U2)
        ax2.plot(solCut[:, 0], solCut[:, 1], 'red', lw = 1)

    if png:
        if debug:
            print('Writing ' + filename_out)
        plt.savefig(filename_out, dpi=300, bbox_inches='tight')
        if debug:
            print('Wrote ' + filename_out)

    plt.show()
