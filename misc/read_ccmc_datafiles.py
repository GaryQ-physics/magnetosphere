import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
from util import urlretrieve, tpad
from units_and_constants import phys
import cxtransform as cx


def time2ccmc_datafile(filename):
    return ''


def getdata(filename, debug=True):
    if type(filename) != str:
        filename = time2ccmc_datafile(filename)
    if not os.path.exists(conf['run_path'] + filename):
        urlretrieve(conf['run_url'] + filename, conf['run_path'] + filename)
        if debug:
            print('Downloading' + filename)
    fname = conf['run_path'] + filename

    data = np.genfromtxt(fname, skip_header=6, comments=None)
    headers = np.loadtxt(fname, dtype=str, skiprows=4, max_rows=1, comments=None)
    assert(headers[0] == '#')
    headers = headers[1:]
    if debug:
        print(headers.shape)
        print(data.shape)

    return [data, headers]

def writevtk(filename, debug=True):
    data, headers = getdata(filename, debug=debug)

    print(headers[7]+headers[8]+headers[9])
    X = data[:, 7:10]*phys['m']
    print(X.shape)

    fname = 'read_ccmc_datafile_grid.vtk'

    print("Writing " + fname)
    f = open(fname,'w')

    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk output\n')
    #f.write('ASCII\n')
    f.write('BINARY\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('\n')
    f.write('POINTS ' + str(X.shape[0]) + ' float\n')
    #np.savetxt(f, X)
    X = np.array(X, dtype='>f')
    f.write(X.tobytes())
    f.close()
    print("Wrote " + fname)

def MAGtoSM(v_MAG, time, ctype_in, ctype_out):
    return cx.transform(v_MAG, time, 'MAG', 'SM', ctype_in=ctype_in, ctype_out=ctype_out)

def MAGtoGEI(v_MAG, time, ctype_in, ctype_out):
    return cx.transform(v_MAG, time, 'MAG', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out)

def GEOtoGEI(v_GEO, time, ctype_in, ctype_out):
    return cx.transform(v_GEO, time, 'GEO', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out)

def GEOtoSM(v_GEO, time, ctype_in, ctype_out):
    return cx.transform(v_GEO, time, 'GEO', 'SM', ctype_in=ctype_in, ctype_out=ctype_out)

def explore(ccmcSite=True):
    data, headers = getdata('YKC_pointdata_754297001306.txt')
    time = np.array(data[:, 0:7], dtype=int)
    #print(headers)
    X = data[:, 7:10]*phys['m']
    if ccmcSite:
        MLAT = 68.8576
        MLON = 300.063
    else:
        MLAT=68.93
        MLON = 299.36

    station_MAG = np.array([1., MLAT, MLON])

    X_conv = MAGtoSM(station_MAG, time, 'sph', 'car')
    print(X_conv)
    print(X)
    #R = np.sqrt(X_conv[:,0]**2 + X_conv[:,1]**2 + X_conv[:,2]**2)
    R = np.sqrt(np.einsum('ij,ij->i', X, X))
    print(R)
    print( np.abs(R-1.) <= 1e-4 )

    assert(np.all(np.abs(R-1.) <= 1e-4))
    assert(np.all(np.abs(X_conv[:,0]-X[:,0]) <= 1e-2))
    assert(np.all(np.abs(X_conv[:,1]-X[:,1]) <= 1e-2))
    assert(np.all(np.abs(X_conv[:,2]-X[:,2]) <= 1e-2))

    import matplotlib.pyplot as plt
    if False:
        plt.plot(X_conv[:,0], 'b.', label='x_conv')
        plt.plot(X_conv[:,1], 'r.', label='y_conv')
        plt.plot(X_conv[:,2], 'g.', label='z_conv')
    if True:
        plt.plot(X_conv[:,0]-X[:,0], 'b.', label='xdif')
        plt.plot(X_conv[:,1]-X[:,1], 'r.', label='ydif')
        plt.plot(X_conv[:,2]-X[:,2], 'g.', label='zdif')
    if False:
        plt.plot(X[:,0], 'b.', label='X')
        plt.plot(X[:,1], 'r.', label='Y')
        plt.plot(X[:,2], 'g.', label='Z')
    if False:
        plt.plot(R-1., 'k.', label='Rdif')
    
    plt.xlabel('index')
    plt.legend()

    #plt.ylabel('X')
    #plt.axvline(x=1.25)
    plt.grid()
    plt.show()


def plot(fileyear, pltvars=['dBn', 'facdBn', 'sumBn', 'JhdBn', 'JpBn']):
    import matplotlib.pyplot as plt
    #import matplotlib.dates
    import datetime

    data, headers = getdata(str(fileyear) + '_YKC_pointdata.txt')
    time = np.array(data[:, 0:7], dtype=int)
    assert(np.all(time[:,0] == fileyear))

    dtimes = []
    for i in range(time.shape[0]):
        dtimes.append(datetime.datetime(time[i,0],time[i,1],time[i,2],time[i,3],time[i,4],time[i,5]))

    #dtimes = matplotlib.dates.date2num(dtimes)

    for var in pltvars:
        toSum = var.split('+')
        values = np.zeros((data.shape[0],))
        for v in toSum:
            ind = list(headers).index(v)
            values = values + data[:,ind]
        plt.plot(dtimes, values, label=var)
    '''
    for i in range(data.shape[1]):
        if headers[j] in pltvars:
            #plt.plot_date(
            plt.plot(dtimes, data[:,j], label=headers[j])
    '''
    # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
    plt.gcf().autofmt_xdate()

    plt.xlabel('year ' + str(fileyear))
    plt.legend()
    plt.show()

plot(2006, pltvars = ['dBn', 'facdBn'])
plot(2006, pltvars = ['JhdBn', 'JpBn', 'JhdBn+JpBn'])
plot(2006, pltvars = ['sumBn'])

"""
working hypothesis:
    'JhdBn+JpBn' = ionosphere contribution to deltaB in fig3a in CalcDeltaB
    'facdBn' = FAC contributions to deltaB in fig2a in CalcDeltaB

    'dBn' is magnetosphere contribution output by SWMF, see 
    https://ccmc.gsfc.nasa.gov/VIS_DOCS/MAGNETOMETER_timeline_quantities.php

"""

"""
---- https://ccmc.gsfc.nasa.gov/VIS_DOCS/Magnetometer_names.php ----
 Station                 Code    MLAT         MLON
 Yellowknife             YKC      68.8576     300.063




---- CalcDeltaB paper ----
  Station                Geographic            Geomagnetic
Name   IAGA_Code    Latitude  Longitude    Latitude  Longitude

Yellowknife YKC      62.48     245.52       68.93     299.36


"""
