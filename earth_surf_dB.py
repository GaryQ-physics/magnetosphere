import os
import sys
import numpy as np
import tempfile


from config import conf
import biot_savart_kameleon_interpolated_grid as bsk
import util
import cxtransform as cx

def plot(LON, LAT, data, title):
    import matplotlib.pyplot as plt

    fig = plt.figure()
    axes = fig.gca()

    #axes.set_title(title, fontsize=10, family='monospace')

    Nbt = 4 # Approximate number of colors between ticks for linear scale
    import matplotlib
    cmap = matplotlib.pyplot.get_cmap('viridis', Nbt*(1000-1))

    axes.set(title=title)
    axes.set(xlabel='geo lon')
    axes.set(ylabel='geo lat')
    axes.axis('square')
    axes.set_xlim(-180., 180.)
    axes.set_ylim(-90., 90.)

    datamax = np.max(data)
    datamin = np.min(data)
    print(datamax,datamin)

    import matplotlib.colors as colors
    norm = colors.SymLogNorm(linthresh=1., vmin=datamin, vmax=datamax)
    pcm = axes.pcolormesh(LON, LAT, data, norm=norm, cmap=cmap, vmin=datamin, vmax=datamax)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    cb = axes.figure.colorbar(pcm, cax=cax)

    plt.show()


def tofile(run, time, para=False):
    xlims=(-16., 16.); ylims=(-16., 16.); zlims=(-16., 16.); d=0.25
    #xlims=(-48., 16.); ylims=(-32., 32.); zlims=(-32., 32.); d=0.25

    Nlat = 19
    Nlon = 37

    R = 1.
    lat = np.linspace(-90., 90., Nlat)
    lon = np.linspace(-180., 180., Nlon)

    LON, LAT = np.meshgrid(lon, lat)
    LAT = LAT.flatten()
    LON = LON.flatten()

    earthgrid = np.column_stack([np.ones(Nlon*Nlat), LAT, LON])

    locations = cx.GEOtoMAG(earthgrid, time, 'sph', 'sph')

    surfB_north = np.zeros(locations.shape[0])
    surfB_east = np.zeros(locations.shape[0])
    surfB_down = np.zeros(locations.shape[0])

    def comp(i):
        deltaB = bsk.integrate(run, time, locations[i,1], locations[i,2], para=False,
            xlims=xlims, ylims=ylims, zlims=zlims, d=d, returnAll=False)
        deltaB_loc = bsk.toMAGLocalComponents(time, locations[i,1], locations[i,2], deltaB)
        print('i = ' + str(i))
        print('deltaB_loc = ' + str(deltaB_loc))
        return deltaB_loc

    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        surfB = Parallel(n_jobs=num_cores)(\
                    delayed(comp)(i) for i in range(Nlat*Nlon))
    else:
        surfB = []
        for i in range(Nlat*Nlon):
            surfB.append(comp(i))


    surfB = np.array(surfB)
    surfB_north = surfB[:,0]
    surfB_east = surfB[:,1]
    surfB_down = surfB[:,2]

    def safenumpy_tofile(array, outname):
        fname = outname
        while os.path.exists(fname):
            fname = fname + '-old.bin'
        os.system('mv %s %s'%(outname,fname))
        array.tofile(outname)

    subdir = '%04d%02d%02dT%02d%02d%02d/' % tuple(util.tpad(time, length=6))
    direct = conf[run+'_derived'] + subdir

    if not os.path.exists(direct):
        os.makedirs(direct)

    print('writing files')
    safenumpy_tofile(surfB_north, direct + 'surfB_north.bin')
    safenumpy_tofile(surfB_east,  direct + 'surfB_east.bin')
    safenumpy_tofile(surfB_down,  direct + 'surfB_down.bin')
    safenumpy_tofile(LON, direct + 'LON.bin')
    safenumpy_tofile(LAT, direct + 'LAT.bin')
    print('wrote files')

    surfB_north = surfB_north.reshape((Nlat, Nlon))
    surfB_east = surfB_east.reshape((Nlat, Nlon))
    surfB_down = surfB_down.reshape((Nlat, Nlon))

    LON = LON.reshape((Nlat, Nlon))
    LAT = LAT.reshape((Nlat, Nlon))

    print(surfB_north)
    print(LAT)
    print(LON)


    #plot(LON, LAT, surfB_north)

def fromfile(run, time):
    Nlat = 19
    Nlon = 37

    subdir = '%04d%02d%02dT%02d%02d%02d/' % tuple(util.tpad(time, length=6))
    direct = conf[run+'_derived'] + subdir

    surfB_north = np.fromfile(direct + 'surfB_north.bin')
    surfB_east = np.fromfile(direct + 'surfB_east.bin')
    surfB_down = np.fromfile(direct + 'surfB_down.bin')
    LON = np.fromfile(direct + 'LON.bin')
    LAT = np.fromfile(direct + 'LAT.bin')

    surfB_north = surfB_north.reshape((Nlat, Nlon))
    surfB_east = surfB_east.reshape((Nlat, Nlon))
    surfB_down = surfB_down.reshape((Nlat, Nlon))

    surfB_Magnitude = np.sqrt(surfB_north**2 + surfB_east**2 + surfB_down**2)

    LON = LON.reshape((Nlat, Nlon))
    LAT = LAT.reshape((Nlat, Nlon))

    print(surfB_north)
    print(LAT)
    print(LON)


    plot(LON, LAT, surfB_east, title='$deltaB_E$')
    plot(LON, LAT, surfB_Magnitude, title='$deltaB_{Magnitude}$')


if __name__ == '__main__':
    #tofile('TESTANALYTIC', (2000,1,1,1,1,0), para=True)
    fromfile()
