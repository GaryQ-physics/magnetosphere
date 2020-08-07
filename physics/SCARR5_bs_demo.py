"""
whern para is false:
    with N = 2560+1 (so dx=dy=dz=0.1) , it takes about 48s to run through each slice.
        Thus its estimated running to completion would take 34.1 hours.
        Assuming the 24 core comp in parallel decreases by factor 1/24, then 1.4hours.
    with N = 64+1 (so dx=dy=dz=4.), it prints Btot_EW = -28.53599693175405
    with N = 128+1 (so dx=dy=dz=2.), it prints Btot_EW = -13.94200518544455

Typical outputs:

    Nx, Ny, Nz = 65, 33, 33
    fullVolume = False
    spacepy_like = False
    slices = 65
    points in slice = 1089
    total points = 70785
    X[0], X[-1], dx = -96.000000, 32.000000, 2.000000
    Y[0], Y[-1], dy = -32.000000, 32.000000, 2.000000
    Z[0], Z[-1], dz = -32.000000, 32.000000, 2.000000
    para = True
    time to process all slices (not including suming up) = 0.13764 min
    Btot = 
    [-2.3255107  -1.41573583 -7.77186095]
    Btot_norm = 8.234933548023848
"""

import numpy as np
import biot_savart_kameleon_interpolated_grid as bsk

global_x_min = -224.
global_x_max = 32.
global_y_min = -128.
global_y_max = 128.
global_z_min = -128.
global_z_max = 128.

data = np.array([[2003, 11, 20, 7, 0, 57.50, 176.00]])
#data = events()

time = data[0, 0:5]
mlat = data[0, 5]
mlon = data[0, 6]
Event = data[0, :]

para = True
fullVolume = False
spacepy_like = False


n = 128 # for testing
#n = 2560

bsk.integrate(time, mlat, mlon, para=para,
    fullVolume=fullVolume, spacepy_like=spacepy_like,
    Nx=n/2+1, xlims=(global_x_min + 128., global_x_max), dx=None,
    Ny=n/4+1, ylims=(global_y_min/4, global_y_max/4), dy=None,
    Nz=n/4+1, zlims=(global_z_min/4, global_z_max/4), dz=None,
    print_output=True)
