# coding=utf-8
import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from units_and_constants import phys
import biot_savart as bs
from biot_savart_demo1 import J_kunits
from probe import probe
from util import time2filename, maketag
import cxtransform as cx



import biot_savart_kameleon_interpolated_grid as bsk

def cdf_to_structured_grid(run, time, mlat, mlon, var, para=True,
        xlims=(-56., 8.), ylims=(-32., 32.), zlims=(-32., 32.),
        d=0.125):

    if dB in var:

         = bsk.run(time, mlat, mlon, filename = conf[run + '_cdf'] + time2filename_ext(run, time),
            para=para, xlims=xlims, ylims=ylims, zlims=zlims, d=d,
            print_output=True, tonpfile=True):

    else:
        x0 = cx.MAGtoGSM(np.array([1., mlat, mlon]), time, 'sph', 'car')

        ret = bs.make_grid(xlims, ylims, zlims, d, d, d)
        Xgrid = ret[0]

        Aa = probe(conf[run + '_cdf'] + time2filename_ext(run, time), Xgrid, var=var, debug=debug)

    out_filename = ' bla bla ' #!!!!!!!!!!!
    writevtk(out_filename, Xgrid, Aa, [Nx,Ny,Nz], 'SCALARS',
             point_data_name = var, title='structured grid ' + run, ftype='BINARY', grid='STRUCTURED_GRID', debug=True):





rbody = 1.25 #???
global_x_min = -224.
global_x_max = 32.
global_y_min = -128.
global_y_max = 128.
global_z_min = -128.
global_z_max = 128.
# difference of 256 for all

# using global min max and dx=dy=dz=0.1  -->  1.6777216e+10 grid points (16 billion)

# 1024 x 1024 x 1024 grid -> 1.073741824 billion  np.zeros works  (i.e.  +/-32 at 0.0625 or +/-64 at 0.125 or  +/-128(full) at 0.25)
# 2048 x 2048 x 2048 grid -> 8.589934592 billion  np.zeros doesnt work


Test = False
debug = False

r_min = 3. # this is minimum distance for kameleon built in tracer to work
           # and also one of nums cited in CalcDeltaB (Lutz Rastatter, Gabor Toth, Maria M. Kuznetsova, Antti A. Pulkkinen)

def Compute(Event, var, calcTotal=False, retTotal=False, dx=.3, dy=.3, dz=.3):
    """
    Either:
        Event = [year, month, day, hours, minutes, MLAT, MLON]    (MLAT and MLON in degrees)
    or:
        Event = [[r, MLAT, MLON], time_list]

    Given an event time (and location), outputs an array for a grid of GSM coordinates
        covering the magnetosphere, and a corresponding array of the physical quantity
        var at those grid points.

    Note: the location part of Event is only used if var is a quantity measured
        relative to the event location, for example dB.

    var can be one of the following strings: 'J', 'dB_EW', 'dB_NS', 'bx', 'by', 'bz', 'jx', 'jy', 'jz', 'p', 'rho'
    """

    if isinstance(Event[0],list):
        mag = np.array(Event[0])
        time = Event[1]
    else:
        time = Event[0:5]
        mag = np.array([1., Event[5], Event[6]])

    if Test:
        X0 = np.array([0., 0.75, 0.])
        Npole = np.array([0., 0., 1.])
    else:
        X0 = cx.MAGtoGSM(mag, time, 'sph', 'car')
        Npole = cx.GEOtoGSM([0., 0., 1.], time, 'car', 'car')

    # datafile
    filename = time2filename(time)
    #conf["run_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-%03d' % tuple(time) + '.out.cdf'

    #X = np.concatenate((np.arange(-200,-20.05,dx_tail), np.arange(-20.,15.,dx) ))
    #tail = -200.
    #tail = -100.
    tail = -20.

    ret = bs.make_grid([tail, 15.], [-30., 30.], [-30., 30.], dx, dy, dz)
    Xgrid = ret[0]
    
    print('X0=',X0)
    if 'dB' in var or var=='J':
        unit_v = np.array([1., 0., 0.])
        if '_' in var:
            a2 = np.cross(Npole, X0)
            a1 = np.cross(X0, a2)
            a1 = a1/np.linalg.norm(a1)
            a2 = a2/np.linalg.norm(a2)
            if '_EW' in var:
                unit_v = a2
            if '_NS' in var:
                unit_v = a1                
        unit_v = np.repeat([unit_v], Xgrid.shape[0], axis=0)

        if Test:
            Jin = J_kunits(Xgrid)*(phys['muA']/phys['m']**2)
        else:
            Jin = probe(time, Xgrid, ['jx','jy','jz'])*(phys['muA']/phys['m']**2)
        if debug:
            print('unit_v=',unit_v)
            print('Jin=',Jin)
            print('Jin sum=', np.sum(Jin))

        if var=='J':
            Aa = Jin
        else:
            dB = bs.deltaB('dB', X0, Xgrid, Jin, V_char = dx*dy*dz) #/nT
            if var == 'dB':
                Aa = np.sqrt(np.sum(dB**2, axis=1))
            else:
                # https://stackoverflow.com/questions/15616742/vectorized-way-of-calculating-row-wise-dot-product-two-matrices-with-scipy
                Aa = np.einsum('ij,ij->i', dB, unit_v) 
        if debug:
            print('Aa=',Aa)
    else:
        Aa = probe(time, Xgrid, var=var, debug=debug)

    if calcTotal:
        total = np.sum(Aa)
        print('Test=',Test)
        print('dx,dy,dz = ', dx, dy, dz)
        print('total = ', total)
        print(Xgrid.shape)
        print(Aa.shape)
        if retTotal:
            return total
    return [Aa, Xgrid, ret[1], ret[2], ret[3]]

def writevtk(out_filename, points, point_data, connectivity, texture,
             point_data_name = 'point data', title='Title', ftype='BINARY', grid='STRUCTURED_GRID', debug=True):
    """
    ftype = 'BINARY' or 'ASCII'
    grid = 'STRUCTURED_GRID' or 'TRIANGLE_STRIPS' or
    """
    if grid == 'STRUCTURED_GRID':
        Nx,Ny,Nz = connectivity

    assert(out_filename[0] == '/')
    f = open(out_filename,'w')
    if debug:
        print("Writing " + fname)
    f.write('# vtk DataFile Version 3.0\n')
    f.write(title + '\n')
    f.write(ftype + '\n')
    f.write('DATASET ' + grid + '\n')

    if grid == 'STRUCTURED_GRID':
        f.write('DIMENSIONS ' + str(Nx) + ' ' + str(Ny) + ' ' + str(Nz) + '\n' )
        f.write('POINTS '+str(Nx*Ny*Nz)+' float\n')

    if ftype=='BINARY':
        Bb = np.array(Bb, dtype='>f')
        f.write(Bb.tobytes())
    elif ftype=='ASCII':
        np.savetxt(f, Bb)

    f.write('\n')
    if grid == 'STRUCTURED_GRID':
        f.write('POINT_DATA ' + str(Nx*Ny*Nz) + '\n')

    f.write(texture + ' ' + point_data_name + ' float\n') # number with float???
    # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping

    '''
    f.write('TEXTURE_COORDINATES TextureCoordinates 2 float\n') # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping


    if var=='J':
        f.write('VECTORS ' + var + ' float\n')
    else:
        f.write('SCALARS ' + var + ' float 1\n')
        f.write('LOOKUP_TABLE default\n')
    '''

    if ftype=='BINARY':
        Aa = np.array(Aa, dtype='>f')
        f.write(Aa.tobytes())
    elif ftype=='ASCII':
        np.savetxt(f, Aa)

    f.close()
    if debug:
        print("Wrote " + fname)

'''

def writevtk(Event, var, calcTotal=False, binary=True, dx=.3, dy=.3, dz=.3, fname=None):
    Aa, Bb, Nx, Ny, Nz = Compute(Event, var, calcTotal=calcTotal, dx=dx, dy=dy, dz=dz)

    if isinstance(Event[0],list):
        time = Event[1]
    else:
        time = Event[0:5]

    #tag = '_%04d:%02d:%02dT%02d:%02d:%02d.%03d' % tuple(time)
    tag = maketag(time)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])

    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)

    if fname == None:
        fname = conf["run_path_derived"] + subdir + 'structured_grid_' + var + tag + '.vtk'


    if debug:
        print('also Nx = ',Nx)

    f = open(fname,'w')
    print("Writing " + fname)
    f.write('# vtk DataFile Version 3.0\n')
    f.write('Structured Grid ' + var + '\n')
    if binary:
        f.write('BINARY\n')
    else:
        f.write('ASCII\n')
    f.write('DATASET STRUCTURED_GRID\n')
    f.write('DIMENSIONS ' + str(Nx) + ' ' + str(Ny) + ' ' + str(Nz) + '\n' )
    f.write('POINTS '+str(Nx*Ny*Nz)+' float\n')
    if binary:
        Bb = np.array(Bb, dtype='>f')
        f.write(Bb.tobytes())
    else:
        np.savetxt(f, Bb)

    f.write('\n')
    f.write('POINT_DATA ' + str(Nx*Ny*Nz) + '\n')
    if var=='J':
        f.write('VECTORS ' + var + ' float\n')
    else:
        f.write('SCALARS ' + var + ' float 1\n')
        f.write('LOOKUP_TABLE default\n')
    if binary:
        Aa = np.array(Aa, dtype='>f')
        f.write(Aa.tobytes())
    else:
        np.savetxt(f, Aa)

    f.close()
    print("Wrote " + fname)
'''
