"""
typical output:

opening: /home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-174600-259.out
---------
dBmhd = 1011.2855377944559
fractional error = -0.20493277341394997
-----------------------
opening: /home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-174500-000.out
---------
dBmhd = 1004.0975310504954
fractional error = -0.24662938647340985


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
import events
import util

import spacepy.pybats.bats as bats

from util import time2filename, dlfile, tpad, urlretrieve

test_Event = np.array([2003, 11, 20, 7, 0, 57.50, 176.00])
#data = np.array([test_Event])
data = events.eventlist()

files = util.filelist()

print files[0]
n = len(files)
n = 2
for i in range(n):
    filename_split = files[i]

    time = util.filename2time(filename_split) # returns length 7 list

    # get list of events with same y,m,d,hour,min
    events_array = events.events(time[0:5])
    print(events_array.shape)
    if events_array.size == 0:
        print('hello there')
        continue

    filename = conf['run_path'] + filename_split[0:34] # remove '.cdf'
    if not os.path.exists(filename):
        print('hello there')
        mess = dlfile(filename, debug=True)
        #urlretrieve(conf['run_url'] + filename_split, filename)

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


    Slice = x==x

    Tr = np.all([-D<=x, x<=D, -D<=y, y<=D, -D<=z, z<=D, Slice], axis=0)
    x_ = x[Tr]
    y_ = y[Tr]
    z_ = z[Tr]
    jx_ = jx[Tr]
    jy_ = jy[Tr]
    jz_ = jz[Tr]

    J_sp = np.column_stack([jx_, jy_, jz_])
    X = np.column_stack([x_, y_, z_])

    for j in range(events_array.shape[0]):
        mlat = events_array[j, 5]
        mlon = events_array[j, 6]
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

'''
    if tuple(Event) == tuple(test_Event):
        from probe import probe
        J_kam = probe(time, X, var=['jx', 'jy', 'jz'], usekV=True)

        out_kam1 = bs.deltaB('deltaB', x0, X, J_kam*(phys['muA']/phys['m']**2), V_char = dV)

        Grid = bs.make_grid(xax, yax, zax, 0, 0, 0)[0]
        J_kam_grid = probe(time, Grid, var=['jx', 'jy', 'jz'], usekV=True)
        out_kam2 = bs.deltaB('deltaB', x0, Grid, J_kam_grid*(phys['muA']/phys['m']**2), V_char = dV)

        L = 4.
        custGrid = bs.make_grid([-L, L], [-L, L], [-L, L], 0.1, 0.1, 0.1)[0]
        J_kam_cust = probe(time, custGrid, var=['jx', 'jy', 'jz'], usekV=True)
        out_kam3 = bs.deltaB('deltaB', x0, custGrid, J_kam_cust*(phys['muA']/phys['m']**2), V_char = 0.1**3)

        print(out_spacepy)
        print(np.linalg.norm(out_spacepy))
        print(out_kam1)
        print(np.linalg.norm(out_kam1))
        print(out_kam2)
        print(np.linalg.norm(out_kam2))
        print(out_kam3)
        print(np.linalg.norm(out_kam3))
'''
