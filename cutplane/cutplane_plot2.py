import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from cutplane2 import data_in_plane, unitvector
from fieldlines import fieldlines
import util
from make_grid import make_grid

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


def plot(run, time, parameter, plane,
         field_lines=None,
         axes=None,
         logz=False, zlims=None, 
         xticks=None, yticks=None, zticks=None,
         dpi=100, showplot=True,
         png=True, pngfile=None, debug=False, **kwargs):#mlat_dB=None, mlon_dB=None):

    """
    # Plot p in x-z plane
    plot(time, 'p', 'xz') 

    # Plot p in plane determined by field line starting at a point in MAG
    # coordinate system; mag = [r, lat, long]
    plot(time, 'p', mag)       

    # Plot p in plane_vs plane
    plot(time, 'p', plane_vs) 
    """
    if isinstance(plane, str):
        plane_str = plane
        plane = {}
        if   plane_str == 'xy':
            xlabel = '$X_{GSM}$ [$R_E$]'
            ylabel = '$Y_{GSM}$ [$R_E$]'
            U1 = np.array([1, 0, 0])
            U2 = np.array([0, 1, 0])
        elif plane_str == 'xz':
            xlabel = '$X_{GSM}$ [$R_E$]'
            ylabel = '$Z_{GSM}$ [$R_E$]'
            U1 = np.array([1, 0, 0])
            U2 = np.array([0, 0, 1])
        elif plane_str == 'yz':
            xlabel = '$Y_{GSM}$ [$R_E$]'
            ylabel = '$Z_{GSM}$ [$R_E$]'
            U1 = np.array([0, 1, 0])
            U2 = np.array([0, 0, 1])
        else:
            raise ValueError ('not a valid coordinate plane')
        plane['e1'] = U1
        plane['e2'] = U2
        plane['e3'] = np.cross(U1, U2)
        plane['ulims'] = (-10,10)
        plane['vlims'] = (-10,10)
        plane['du'] = 0.1
        plane['dv'] = 0.1
    else:
        plane_str = None
        if not isinstance(plane, dict): raise ValueError ('plane must be dict or str')

    try:
        xlabel = kwargs['xlabel']
        ylabel = kwargs['ylabel']
    except KeyError:
        if plane_str is None:
            xlabel = 'u_axis [$R_E$]'
            ylabel = 'v_axis [$R_E$]'

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

    title = run + ' %04d-%02d-%02dT%02d:%02d:%02d.%03d' % tuple(time)
    if 'dB' in parameter:
        assert('mlat_dB' in kwargs.keys())
        assert('mlon_dB' in kwargs.keys())
        import cxtransform as cx
        title = title + '\nat mlat=%.3f, mlon=%.3f, MLT=%.3f hrs'%(kwargs['mlat_dB'], kwargs['mlon_dB'], cx.MAGtoMLT(kwargs['mlon_dB'], time))
 
    filename = util.time2CDFfilename(run, time)

    if pngfile is None:
        outdir = conf[run + '_derived'] + "cutplanes/"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        filename_out = filename.replace(conf[run + 'cdf'], outdir) + '.png'
    else:
        filename_out = pngfile

    x_1d = np.arange(plane['ulims'][0], plane['ulims'][1] + plane['du'], plane['du'])
    y_1d = np.arange(plane['vlims'][0], plane['vlims'][1] + plane['dv'], plane['dv'])
    X, Y = np.meshgrid(x_1d, y_1d) # grid of points on the cut plane
    # note: np.column_stack[X.flatten(), Y.flatten()] == makegrid([x_1d, y_1d, None])[:,:2]
    assert(np.all( np.column_stack([X.flatten(), Y.flatten()]) == make_grid([x_1d, y_1d, None])[:,:2] ))

    if plane_str is None: plane_str = 'na'

    ext = "-" + parameter + '_plane_' + plane_str \
            + '-du_' + str(plane['du']) \
            + '-dv_' + str(plane['dv']) \
            + '-ulims' + str(plane['ulims'][0]) + ',' + str(plane['ulims'][1]) \
            + '-vlims' + str(plane['vlims'][0]) + ',' + str(plane['vlims'][1]) \
            + '.npy'
    cachedir = conf[run + '_derived'] + "cutplanes/cache/" + parameter + "/"
    if 'dB' in parameter:
        cachedir = cachedir + '_%.3f_%.3f/'%(kwargs['mlat_dB'], kwargs['mlon_dB'])
    if not os.path.exists(cachedir):
        os.makedirs(cachedir)
    npfile = filename.replace(conf[run + '_cdf'], cachedir) + ext
    if os.path.exists(npfile):
        print("Reading " + npfile)
        Z = np.load(npfile)
        print("Read " + npfile)
    else:
        if debug:
            print("Interpolating {0:s} onto {1:d}x{2:d} grid" \
                  .format(parameter,len(x_1d),len(y_1d)))
        Z = data_in_plane(run, time, parameter, plane,
                          **kwargs).reshape(X.shape)
        print("Writing " + npfile)
        np.save(npfile, Z)
        print("Wrote " + npfile)

    if field_lines is not None:
        raise RuntimeWarning ('WARNING: need to revisit')
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
        lines = fieldlines(run, time, field_lines)
        for line in lines:
            lineU = np.empty(line.shape)
            for k in range(line.shape[0]):  # TO DO: vectorize this loop
                lineU[k, 0] = np.dot(line[k, :], U1)
                lineU[k, 1] = np.dot(line[k, :], U2)
                lineU[k, 2] = np.dot(line[k, :], U3)
            linesU.append(lineU)

    if run == 'TESTANALYTIC':
            meta = {}
            meta["parameters"] = {
                'jx' : {'plot_unit':'$\\frac{\\mu A}{m^2}$', 'plot_name':'jx'},
                'jy' : {'plot_unit':'$\\frac{\\mu A}{m^2}$', 'plot_name':'jy'},
                'jz' : {'plot_unit':'$\\frac{\\mu A}{m^2}$', 'plot_name':'jz'}
                                    }

    else:
        meta = util.filemeta(filename)

    if 'dB' in parameter:
        #dB_param_dict = {'dB_Magnitude' : '$|\\frac{dB}{dV}|$',
        #                 'dB_north' : '$\\frac{dB_N}{dV}$',
        #                 'dB_east' : '$\\frac{dB_E}{dV}$', 
        #                 'dB_down' : '$\\frac{dB_D}{dV}$'}
        #parameter_unit = '$\\frac{nT}{R_E^3}$'
        dB_param_dict = {'dB_Magnitude' : '$|dB|$',
                         'dB_north' : '$dB_N$',
                         'dB_east' : '$dB_E$', 
                         'dB_down' : '$dB_D$'}
        parameter_unit = '$nT$'
        parameter_label = dB_param_dict[parameter] # get dictionary
    else:
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
    axes.set_xlim(*plane['ulims'])
    axes.set_ylim(*plane['vlims'])

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
