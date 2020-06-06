# cut_plane_python_plotting

import os
import sys
import numpy as np
import pickle
import matplotlib
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
         xlims=[-10,10], ylims=[-10,10], zlims=None, 
         xticks=None, yticks=None, zticks=None,
         dx=0.1, dy=0.1,
         dpi=300, showplot=True,
         png=True, pngfile=None, debug=False):

    """
    plot(time, [r, mlat, mlong], None, 'p')
    plot(time, [GSMx,GSMy,GSMz], [v1, v2], 'p')
    """
    
    xlabel = ''
    ylabel = ''
    if type(pos) == str:
        if pos == 'xy':
            xlabel = 'X [$R_E$]'
            ylabel = 'Y [$R_E$]'
            plane_vs = [[1, 0, 0], [0, 1, 0]]
        if pos == 'xz':
            xlabel = 'X [$R_E$]'
            ylabel = 'Z [$R_E$]'
            plane_vs = [[1, 0, 0], [0, 0, 1]]
        if pos == 'yz':
            xlabel = 'Y [$R_E$]'
            ylabel = 'Z [$R_E$]'
            plane_vs = [[0, 1, 0], [0, 0, 1]]
        pos = [0, 0, 0]
            
    # Plot title
    if plane_vs == None:
        title = 'SCARR5 ' + '%04d%02d%02dT%02d%02d%02d-%03d' % tuple(time)
        title = title + "\n" + "[mlat, mlon]=[{0:.1f}$^o$, {1:.1f}$^o$]".format(pos[1], pos[2])
    else:
        title = 'SCARR5 ' + '%04d-%02d-%02dT%02d:%02d:%02d.%03d' % tuple(time)

    filename = '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-%03d' % tuple(time)

    filename_in = conf["run_path"] + filename + '.out.cdf'
    if pngfile is None:
        filename_out = conf["run_path_derived"] + "cutplanes/" + filename + '.png'
    else:    
        filename_out = pngfile
        
    # open kameleon
    kameleon = ccmc.Kameleon()
    if debug:
        print("Opening " + filename_in)
    kameleon.open(filename_in)
    if debug:
        print("Opened " + filename_in)
    interpolator = kameleon.createNewInterpolator()

    parameter_unit = kameleon.getVisUnit(parameter)
    parameter_unit = parameter_unit.replace('mu', '$\mu$')

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

    x_1d = np.arange(xlims[0], xlims[1] + dx, dx)
    y_1d = np.arange(ylims[0], ylims[1] + dy, dy)
    X, Y = np.meshgrid(x_1d, y_1d) # grid of points on the cut plane
    Z = np.zeros(X.shape)

    if debug:
        print("Interpolating {0:s} onto {1:d}x{2:d} grid".format(parameter,len(x_1d),len(y_1d)))
    #sys.stdout.flush()
    for i in range(X.shape[0]): # note this is y_1d.size (NOT x)
        for j in range(X.shape[1]): 
            # grid of the corresponding values of variable. To be color plotted
            Z[i, j] = data_in_U(kameleon, interpolator, parameter,
                                X[i, j], Y[i, j], U1, U2, U3)

    kameleon.close()
    #sys.stdout.flush()
    if debug:
        print("Closed " + filename_in + "\n")

    r = float(xlims[1]-xlims[0])/float(ylims[1]-ylims[0])
    
    fig = plt.figure(figsize=(6, 6/r), dpi=dpi, tight_layout=False)
    ax2 = fig.gca()

    # Plot cut plane data
    ax2.set_title(title, fontsize=10)
    if plane_vs == None:
        xlabel = "Tailward distance [$R_E$]"
        ylabel = "Northward distance [$R_E$]"
    elif xlabel == '':    
        xlabel = '[%.1f, %.1f, %.1f]' % tuple(U1) + ' dir' + " [$R_E$]"
        ylabel = '[%.1f, %.1f, %.1f]' % tuple(U2) + ' dir' + " [$R_E$]"

    from niceticks import niceticks

    if zticks is not None and zlims is not None:
        zlims = None
        print('Ignoring zlims b/c zticks given.')

    if zlims is not None:
        boundaries = niceticks(zlims[0], zlims[1], 10)
    elif zticks is not None:
        boundaries = zticks
    else:
        boundaries = niceticks(np.min(Z.flatten()), np.max(Z.flatten()), 10)

    if zticks is None:
        zticks = boundaries

    Nbt = 4 # Number of colors between ticks
    cmap = matplotlib.pyplot.get_cmap('viridis', Nbt*len(boundaries))

    ax2.set(xlabel=xlabel)
    ax2.set(ylabel=ylabel)
    ax2.axis('square')
    pcm = ax2.pcolormesh(X, Y, Z, cmap=cmap)

    plt.xlim(xlims[0], xlims[1])
    plt.ylim(ylims[0], ylims[1])

    db = boundaries[1]-boundaries[0]
    boundaries = np.arange(boundaries[0], boundaries[-1] + db/Nbt, db/Nbt)

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax2)
    cax1 = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar(pcm, cax=cax1, ticks=zticks, boundaries=boundaries)
    cb.ax.set_title(parameter + ' [' + parameter_unit + ']', va='bottom', fontsize=10)
    fig.tight_layout(h_pad=0)

    # Add field line to 2D plot
    if plane_vs == None:
        solCut = np.zeros((sol.shape[0], 2))
        for k in range(sol.shape[0]):
            solCut[k, 0] = np.dot(sol[k, :], U1)
            solCut[k, 1] = np.dot(sol[k, :], U2)
        ax2.plot(solCut[:, 0], solCut[:, 1], 'red', lw = 1)

    if png:
        if debug:
            print('Writing ' + filename_out)
        plt.savefig(filename_out, dpi=dpi, bbox_inches='tight')
        if debug:
            print('Wrote ' + filename_out)

    if showplot:
        plt.show()

    info = {'min': np.min(Z.flatten()),
            'max': np.max(Z.flatten())
            }
    
    return info