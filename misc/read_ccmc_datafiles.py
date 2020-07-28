"""
This script reads the calculated values for the different contributions to the
(probably) magnetic field that are in the .txt files from
https://ccmc.gsfc.nasa.gov/results/viewrun.php?domain=GM&runnumber=SWPC_SWMF_052811_2


"""

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
from util import urlretrieve, tpad
from units_and_constants import phys
import cxtransform as cx


def year2ccmc_datafile(tofilename):
    """
    tofilename = [2006, 'YKC']
    """
    return str(tofilename[0]) + '_' + tofilename[1] + '_pointdata.txt'


def getdata(filename, debug=False):
    """

    Parameters
    ----------
    filename : string
            name of .txt datafile (NOT full path)
    debug : boolean, OPTIONAL
        DESCRIPTION. The default is True.

    Returns
    -------
    list of np arrays = [data, headers]
    
    where data is (N,25) array of N datapoints given in the file and 
    headers is (25,) array of string for the corresponding quantities of data

    """

    if type(filename) != str:
        filename = year2ccmc_datafile(filename)


    fname = conf['data_path'] + filename
    data = np.genfromtxt(fname, skip_header=6, comments=None)
    headers = np.loadtxt(fname, dtype=str, skiprows=4, max_rows=1, comments=None)
    assert(headers[0] == '#')
    headers = headers[1:]
    if debug:
        print(data.shape)
        print(headers.shape)

    return [data, headers]

def writevtk(filename, debug=True):
    """

    Parameters
    ----------
    filename : string
            name of .txt datafile (NOT full path)
    debug : boolean, OPTIONAL
        DESCRIPTION. The default is True.

    Returns
    -------
    None.
    
    writes binary UNSTRUCTURED_GRID vtk file of all the cartesian x,y,z points
        given by data in getdata()

    """
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


def explore(fileyear=2006, useGEO=True, magSource=0):
    """
    working hypothesis:
        'JhdBn+JpBn' = ionosphere contribution to deltaB in fig3a in CalcDeltaB
        'facdBn' = FAC contributions to deltaB in fig2a in CalcDeltaB
    
        'dBn' is magnetosphere contribution output by SWMF, see 
        https://ccmc.gsfc.nasa.gov/VIS_DOCS/MAGNETOMETER_timeline_quantities.php
    
        X,Y,Z columns are the station position in SM coordinates at that time
    
    ##########################################################################
    
    ---- https://ccmc.gsfc.nasa.gov/VIS_DOCS/Magnetometer_names.php ----
     Station                 Code    MLAT         MLON
     Yellowknife             YKC      68.8576     300.063
    
    ---- CalcDeltaB paper ----
      Station                Geographic            Geomagnetic
    Name   IAGA_Code    Latitude  Longitude    Latitude  Longitude
    
    Yellowknife YKC      62.48     245.52       68.93     299.36
    
    ---- https://www.geomag.nrcan.gc.ca/obs/ykc-en.php ----
    IAGA alphabetic code 	YKC
    IAGA numeric code 	028246
    Geographic coordinates 	62.480 N, 245.518 E
    Geomagnetic coordinates (IGRF-12 (2015)) 	68.62 N, 58.22 W (2015.0)
    Elevation 	198 m

    ---- http://www.wdc.bgs.ac.uk/obsinfo/ykc.html ----
    Latitude 	62.48AN
    Longitude 	245.518AE
    Altitude 	198.0m

    """
    data, headers = getdata(fileyear)
    time = np.array(data[:, 0:7], dtype=int)
    #print(headers)
    X = data[:, 7:10]*phys['m']

    if useGEO:

        lat = 62.480
        lon = 245.518
        station_GEO = np.array([1., lat, lon])
        X_conv = cx.GEOtoSM(station_GEO, time, 'sph', 'car')

    else:

        if magSource==0:
            MLAT = 68.8576
            MLON = 300.063
        elif magSource==1:
            MLAT=68.93
            MLON = 299.36
        station_MAG = np.array([1., MLAT, MLON])
        X_conv = cx.MAGtoSM(station_MAG, time, 'sph', 'car')


    di_location = cx.dipole(time)
    diLAT = di_location[:,1]
    diLON = di_location[:,2]

    #print(X_conv)
    #print(X)
    #R = np.sqrt(X_conv[:,0]**2 + X_conv[:,1]**2 + X_conv[:,2]**2)
    R = np.sqrt(np.einsum('ij,ij->i', X, X))
    #print(R)
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
    if False:
        plt.plot(X_conv[:,0]-X[:,0], 'b.', label='xdif')
        plt.plot(X_conv[:,1]-X[:,1], 'r.', label='ydif')
        plt.plot(X_conv[:,2]-X[:,2], 'g.', label='zdif')
    if False:
        plt.plot(X[:,0], 'b.', label='X')
        plt.plot(X[:,1], 'r.', label='Y')
        plt.plot(X[:,2], 'g.', label='Z')
    if False:
        plt.plot(R-1., 'k.', label='Rdif')
    if True:
        plt.plot(diLAT, 'b.', label='diLAT')
        plt.plot(diLON, 'b.', label='diLON')
        plt.axvline(x=79.0)
        plt.axvline(x=289.1)

    
    plt.xlabel('index')
    plt.legend()

    #plt.ylabel('X')
    #plt.axvline(x=1.25)
    plt.grid()
    plt.show()
    return data


def plot(fileyear, pltvars=['dBn', 'facdBn', 'sumBn', 'JhdBn', 'JpBn', 'JhdBn+JpBn']):
    """

    Parameters
    ----------
    fileyear : integer or string
        the year the datafile is from, or just the string name of the datafilefile. It is implicitly assumed there is only one 
        datafile stored for a given year.
        so far tried with 2001 and 2006
    pltvars : list, OPTIONAL
        A list of strings of all the variable to be plotted. The list elements
        can either be a variable in the headers of the file or a sum of them.
        when denoting a sum, use no spaces in the string, just +
        
        The default is ['dBn', 'facdBn', 'sumBn', 'JhdBn', 'JpBn', 'JhdBn+JpBn'].

    Returns
    -------
    None.
    
    plots (using matplotlib) a graph of all the variables in pltvars as a function
    of time in datetime format

    """
    import matplotlib.pyplot as plt
    #import matplotlib.dates
    import datetime

    data, headers = getdata(fileyear)
    time = np.array(data[:, 0:7], dtype=int)

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

#plot(2006, pltvars = ['dBn', 'facdBn'])
#plot(2006, pltvars = ['JhdBn', 'JpBn', 'JhdBn+JpBn'])
#plot(2006, pltvars = ['sumBn'])
#explore()
