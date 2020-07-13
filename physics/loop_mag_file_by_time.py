"""
this script loops over different locations on earth surface given a time
however these locations are not nessesarely events in LOCALIZED.txt

when run in terminal, this script prints out series of:
    locations: year, month, day, hour, minuts; latitude in MAG, longitude in MAG
    dBmhd = the magnitude of dBmhd calculated from the components 
            dBnMhd, dBeMhd, dBdMhd as givent in the mag_grid____.out files.
            i.e.  it is the magnitude of dB given in SWMF output due to 
            magnetosphere currents.
    fractional error = the fractional error, relative to dBmhd, of the
                    result gotten by performing biot savart (using the 
                    biot_savart module) on the SWMF fine mesh grid near
                    earth.

typical out:
    see loop_mag_file_by_time.log
"""

import sys
import os
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import numpy as np

from units_and_constants import phys
import cxtransform as cx
import biot_savart as bs
import read_mag_grid_files as rmg
import util
import read_mag_grid_files as rmg

import spacepy.pybats.bats as bats


time = [2003,11,20,7,0]

data, headers = rmg.getdata(time)

n = data.shape[0]
n = 10
for j in range(n):
    i = j + 10000
    mlon = data[i, 0]
    mlat = data[i, 1]
    # note this is oposite convention, but thats what the mag file is

    filename = util.time2filename(time, extension='.out')
    if not os.path.exists(filename):
        print('hello there')
        mess = dlfile(filename, debug=True)

    print('opening: ' + filename)
    # read in the 3d magnetosphere
    data3d = bats.Bats2d(filename)

    # get the cell coordinates
    x = data3d['x']
    y = data3d['y']
    z = data3d['z']
    jx = data3d['jx']
    jy = data3d['jy']
    jz = data3d['jz']

    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    z = np.array(z, dtype=float)
    jx = np.array(jx, dtype=float)
    jy = np.array(jy, dtype=float)
    jz = np.array(jz, dtype=float)

    D = 3.96875
    xax = np.linspace(-D, D, 128)
    yax = np.linspace(-D, D, 128)
    zax = np.linspace(-D, D, 128)
    dx = 0.0625
    dV = dx**3

    assert((D+D)/(128-1) == dx)
    assert(xax[1] - xax[0] == dx)

    Tr = np.all([-D<=x, x<=D, -D<=y, y<=D, -D<=z, z<=D], axis=0)
    x_ = x[Tr]
    y_ = y[Tr]
    z_ = z[Tr]
    jx_ = jx[Tr]
    jy_ = jy[Tr]
    jz_ = jz[Tr]

    J_sp = np.column_stack([jx_, jy_, jz_])
    X = np.column_stack([x_, y_, z_])

    x0 = cx.MAGtoGSM([1., mlat, mlon], time, 'sph', 'car')

    out_spacepy = bs.deltaB('deltaB', x0, X, J_sp*(phys['muA']/phys['m']**2), V_char = dV)
    B_spacepy = np.linalg.norm(out_spacepy)
    magfile = rmg.analyzedata(time, mlat, mlon, debug=False)
    Bmhd_magfile = np.sqrt(magfile[1]**2 + magfile[2]**2 + magfile[3]**2)

    print('---')
    #print(B_spacepy)
    print('Event: {0:d},{1:d},{2:d},{3:d},{4:d};{5:.2f},{6:.2f}'.format\
                (time[0],time[1],time[2],time[3],time[4],mlat,mlon))
    print('dBmhd = ' + str(Bmhd_magfile))
    print('fractional error = ' + str( (B_spacepy-Bmhd_magfile)/Bmhd_magfile ))
    print('----------------------------')
