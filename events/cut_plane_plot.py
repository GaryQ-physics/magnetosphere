# cut_plane_python_plotting

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import _CCMC as ccmc
from cut_plane import ex_data, Compute


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


def plot(time, pos, plane_vs, parameter,
         xlim=[0,4], ylim=[-3,3], dx=0.1, dy=0.1, png=True, debug=False):
    # plot(time, [r, mlat, mlong], None, 'p')
    # plot(time, [GSMx,GSMy,GSMz], [v1, v2], 'p')
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    
    if type(pos) == str:
        if pos == 'xy':
            plot(time, [0, 0, 0], [[1, 0, 0], [0, 1, 0]],
                 parameter, xlim=xlim, ylim=ylim)
            return
        if pos == 'xz':
            plot(time, [0, 0, 0], [[1, 0, 0], [0, 0, 1]],
                 parameter, xlim=xlim, ylim=ylim)
            return
        if pos == 'yz':
            plot(time, [0, 0, 0], [[0, 1, 0], [0, 0, 1]],
                 parameter, xlim=xlim, ylim=ylim)
            return
            
    # Plot title
    if plane_vs == None:
        title = 'SCARR5 ' + '%04d%02d%02dT%02d%02d%02d' % tuple(time)
        title = title + "\n" + "[mlat,mlon]=[{0:.1f}, {1:.1f}]".format(pos[1], pos[2])
    else:
        title = 'SCARR5 ' + '%04d%02d%02dT%02d%02d%02d' % tuple(time)
        title = title + "\n" + "GSM coords"

    filename = '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % tuple(time)

    filename_in = conf["run_path"] + filename + '.out.cdf'
    filename_out = conf["run_path_derived"] + filename + '.png'
    
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

    if plane_vs == None:
        r = pos[0]
        mlat = pos[1]
        mlon = pos[2]
        Event = time + [mlon, mlat]
        if debug: print('Event=',Event)
        ret = Compute(Event, ret_sol=True, r=r)
        U1 = ret[1]
        U2 = ret[2]
        U3 = ret[3]
        sol = ret[4]
        
        if debug: print('sol = ', sol)
    else:
        U1 = np.array(plane_vs[0])
        U2 = np.array(plane_vs[1])
        U3 = np.cross(U1, U2)

    x_1d = np.arange(xlim[0], xlim[1], dx)
    y_1d = np.arange(ylim[0], ylim[1], dy)
    X, Y = np.meshgrid(x_1d, y_1d) # grid of points on the cutplane
    Z = np.zeros(X.shape)
    if debug:
        print x_1d.shape,y_1d.shape
        print X.shape,Y.shape,Z.shape
        print x_1d.size, y_1d.size
        print x_1d
        print y_1d
        print X
    for i in range(X.shape[0]): # note this is y_1d.size (NOT x)
        for j in range(X.shape[1]): 
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
    #xlabel = 'xlable'
    #ylabel = 'ylable'
    #if plane_vs == None and type(pos) == str:
    #    xlabel = pos[0] + '_GSM ' + "[$R_E$]"
    #    ylabel = pos[1] + '_GSM ' + "[$R_E$]"
    if plane_vs == None:
        xlabel = "Tailward distance [$R_E$]"
        ylabel = "Northward distance [$R_E$]"
    else:
        xlabel = '[%.2f, %.2f, %.2f]' %tuple(U1) + ' in GSM' + "[$R_E$]"
        ylabel = '[%.2f, %.2f, %.2f]' %tuple(U2) + ' in GSM' + "[$R_E$]"

    ax2.set(xlabel = xlabel)
    ax2.set(ylabel = ylabel)
    ax2.axis('square')
    pcm = ax2.pcolormesh(X, Y, Z)

    # Reason for choice of fraction and pad:
    # https://stackoverflow.com/a/39948312/1491619
    cb = fig.colorbar(pcm, ax = ax2, fraction = 0.04, pad = 0.08)
    cb.ax.set_title(parameter + ' [$' + parameter_unit + '$]', ha='left')
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])

    # Add field line to 2D plot
    if plane_vs == None:
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
