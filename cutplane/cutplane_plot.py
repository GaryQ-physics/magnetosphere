import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from cutplane import data2d, unitvector, fieldlines
import util


def set_colorbar(ax, pcm, zticks, Nbt, title=None, logz=False):
    import matplotlib.pyplot as plt        
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import  matplotlib.ticker as ticker
    
    def fmt(x, pos):
        # replace("-", u"\u2013") uses a shorter dash for a negative sign
        if abs(x) >= 1e3 or abs(x) <= 1e-3 and x != 0:
            a, b = '{:.0e}'.format(x).split('e')
            b = int(b)
            return r'${} \cdot 10^{{{}}}$'.format(a, b).replace("-", u"\u2013")
        else:
            return r'{0:.1f}'.format(x).replace("-", u"\u2013")

    formatter = ticker.FuncFormatter(fmt)
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2,2))
    
    db = zticks[1] - zticks[0]
    boundaries = np.arange(zticks[0], zticks[-1] + db/Nbt, db/Nbt)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    if logz:
        cb = ax.figure.colorbar(pcm, cax=cax)
    else:         
        cb = ax.figure.colorbar(pcm, cax=cax, ticks=zticks,
                          boundaries=boundaries, format=formatter)
    
    # Force axis labels to be computed and set
    plt.draw()

    # Hide x 10^{n} text which is difficult to position and takes
    # up a lot of space. Put num\cdot 10^{n} on highest tick label.

    ot = cb.ax.yaxis.offsetText.get_text()
    cb.ax.yaxis.offsetText.set_visible(False)
    z = cb.ax.get_yticklabels(0)
    nt = z[-1].get_text()
    nt = nt + ot.replace('times', 'cdot')
    z[-1].set_text(nt)
    cb.ax.set_yticklabels(z)

    if title is not None:
        cb.ax.set_title(title, va='bottom', fontsize=10)

def plot(time, parameter, arg3,
         field_lines=None,
         axes=None,
         logz=False,
         xlims=[-10,10], ylims=[-10,10], zlims=None, 
         xticks=None, yticks=None, zticks=None,
         dx=0.1, dy=0.1,
         dpi=100, showplot=True,
         png=True, pngfile=None, debug=False):

    """
    # Plot p in x-z plane
    plot(time, 'p', 'xz') 

    # Plot p in plane determined by field line starting at a point in MAG
    # coordinate system; mag = [r, lat, long]
    plot(time, 'p', mag)       

    # Plot p in plane_vs plane
    plot(time, 'p', plane_vs) 
    """

    import matplotlib
    if not showplot:
        matplotlib.use("Agg", warn=False)
    if axes is None:
        # This import statement must follow the .use call above.
        import matplotlib.pyplot as plt

    # Set monospace font for changing text
    # TODO: Use with ... to set these.
    matplotlib.rc('font', family='serif')
    matplotlib.rcParams['mathtext.fontset'] = 'dejavuserif'
    
    xlabel = ''
    ylabel = ''
    title = 'SCARR5 ' + '%04d-%02d-%02dT%02d:%02d:%02d.%03d' % tuple(time)
    if type(arg3) == str:
        if arg3 == 'xy':
            xlabel = 'X [$R_E$]'
            ylabel = 'Y [$R_E$]'
            U1 = np.array([1, 0, 0])
            U2 = np.array([0, 1, 0])
        if arg3 == 'xz':
            xlabel = 'X [$R_E$]'
            ylabel = 'Z [$R_E$]'
            U1 = np.array([1, 0, 0])
            U2 = np.array([0, 0, 1])
        if arg3 == 'yz':
            xlabel = 'Y [$R_E$]'
            ylabel = 'Z [$R_E$]'
            U1 = np.array([0, 1, 0])
            U2 = np.array([0, 0, 1])
    else:
        arg3 = np.array(arg3)
        if arg3.size == 3:
            title = title + "\n" + "[mlat, mlon]=[{0:.1f}$^o$, {1:.1f}$^o$]" \
                    .format(arg3[1], arg3[2])
            U = unitvector(time, arg3)
            U1 = U[0]
            U2 = U[1]
        else:
            U1 = arg3[0]
            U2 = arg3[1]
    
        # Normalize the U vectors.
        U1 = U1/np.linalg.norm(U1)
        U2 = U2/np.linalg.norm(U2)
        xlabel = '[%.1f, %.1f, %.1f]' % tuple(U1) + ' dir' + " [$R_E$]"
        ylabel = '[%.1f, %.1f, %.1f]' % tuple(U2) + ' dir' + " [$R_E$]"

    U3 = np.cross(U1, U2)
 
    filename = util.time2filename(time)

    if pngfile is None:
        outdir = conf['run_path_derived'] + "cutplanes/"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        filename_out = filename.replace(conf['run_path'], outdir) + '.png'
    else:
        filename_out = pngfile

    x_1d = np.arange(xlims[0], xlims[1] + dx, dx)
    y_1d = np.arange(ylims[0], ylims[1] + dy, dy)
    X, Y = np.meshgrid(x_1d, y_1d) # grid of points on the cut plane


    if type(arg3) == str:
        ext = "-" + parameter + '_plane_' + arg3 \
                + '-dx_' + str(dx) \
                + '-dy_' + str(dy) \
                + '-xlims' + str(xlims[0]) + ',' + str(xlims[1]) \
                + '-ylims' + str(ylims[0]) + ',' + str(ylims[1]) \
                + '.npy'
        cachedir = conf['run_path_derived'] + "cutplanes/cache/" + parameter + "/"
        if not os.path.exists(cachedir):
            os.makedirs(cachedir)
        npfile = filename.replace(conf['run_path'], cachedir) + ext
        if os.path.exists(npfile):
            print("Reading " + npfile)
            Z = np.load(npfile)
            print("Read " + npfile)
        else:
            if debug:
                print("Interpolating {0:s} onto {1:d}x{2:d} grid" \
                      .format(parameter,len(x_1d),len(y_1d)))
            Z = data2d(time, parameter, X, Y, [U1, U2, U3], debug=debug)
            print("Writing " + npfile)
            np.save(npfile, Z)
            print("Wrote " + npfile)

    if field_lines is not None:
        linesU = []
        '''
        for field_line in field_lines:
            line = fieldlines(time, field_line)
            lineU = np.empty(line.shape)
            for k in range(line.shape[0]):
                lineU[k, 0] = np.dot(line[k, :], U1)
                lineU[k, 1] = np.dot(line[k, :], U2)
                lineU[k, 2] = np.dot(line[k, :], U3)
            linesU.append(lineU)
        '''
        lines = fieldlines(time, field_lines)
        for line in lines:
            lineU = np.empty(line.shape)
            for k in range(line.shape[0]):  # TO DO: vectorize this loop
                lineU[k, 0] = np.dot(line[k, :], U1)
                lineU[k, 1] = np.dot(line[k, :], U2)
                lineU[k, 2] = np.dot(line[k, :], U3)
            linesU.append(lineU)

    meta = util.filemeta(filename)

    parameter_unit = meta["parameters"][parameter]['plot_unit']
    parameter_label = meta["parameters"][parameter]['plot_name']
    
    if axes is None:
        fig = plt.figure(figsize=(7, 6), dpi=dpi, tight_layout=False)
        axes = fig.gca()
        axes_given = False
    else:
        axes_given = True
        if png:
            print("png=True ignored because axes was given as keyword argument.")
        png = False
        showplot = False

    # Plot cut plane data
    axes.set_title(title, fontsize=10, family='monospace')

    from niceticks import niceticks

    if zticks is not None and zlims is not None:
        zlims = None
        print('Ignoring zlims b/c zticks given.')

    if zticks is None and zlims is None:
        zticks = niceticks(np.min(Z.flatten()), np.max(Z.flatten()), 10)
    elif zticks is None and zlims is not None: 
        zticks = niceticks(zlims[0], zlims[1], 10)              

    Nbt = 4 # Approximate number of colors between ticks for linear scale
    cmap = matplotlib.pyplot.get_cmap('viridis', Nbt*(len(zticks)-1))

    axes.set(xlabel=xlabel)
    axes.set(ylabel=ylabel)
    axes.axis('square')
    axes.set_xlim(xlims[0], xlims[1])
    axes.set_ylim(ylims[0], ylims[1])

    if logz:
        import matplotlib.colors as colors
        norm = colors.SymLogNorm(linthresh=0.001, vmin=zticks[0], vmax=zticks[-1])
        pcm = axes.pcolormesh(X, Y, Z, norm=norm, cmap=cmap, vmin=zticks[0], vmax=zticks[-1])
    else:
        pcm = axes.pcolormesh(X, Y, Z, cmap=cmap, vmin=zticks[0], vmax=zticks[-1])

    set_colorbar(axes, pcm, zticks, Nbt,
                 title=parameter_label + ' [' + parameter_unit + ']',
                 logz=logz)
    
    if field_lines is not None:
        for lineU in linesU:
            axes.plot(lineU[:, 0], lineU[:, 1], 'red', lw = 1)

    if axes_given == False:
        if png:
            if debug:
                print('Writing ' + filename_out)
            plt.savefig(filename_out, dpi=dpi) # bbox_inches='tight')
            if debug or pngfile is None:
                print('Wrote ' + filename_out)
    
        if showplot:
            plt.show()
        else:
            plt.close()

    info = {'min': np.min(Z.flatten()),
            'max': np.max(Z.flatten())
            }
    
    return info
