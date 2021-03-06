"""
This script reads the calculated values for the different contributions to the
magnetic field that are in mag_grid____.out files


NOTE:

lats = np.unique(data[:,1])
  -> array([-87.5       , -86.49425287, -85.48850575, -84.48275862,
            ...
            85.48850575,  86.49425287,  87.5       ])
lats[1]-lats[0]
  -> 1.0057471260000028
lats[6]-lats[5]
  -> 1.0057471270000065
lats.shape
  -> (175,)
87.5-(-87.5)
  -> 175.0
175./(175-1)
  -> 1.0057471264367817

lons = np.unique(data[:,0])
  -> array([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
            ...
            352., 353., 354., 355., 356., 357., 358., 359.])
lons.shape
  -> (360,)
(359.-0.)/(360-1)
  -> 1.0

so while lons in steps of 1. , lats are insteps of 175./174 
(to some decimal aprox that is not completely consistent between elements)
so someone probably screwed up the linspace and meant to use (176,)


"""

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
from util import dlfile, tpad


def time2mag_grid_file(time):
    """
    takes in tuple/list/array of time
    returns mag .out filename associated with that time
    note: mag file names only have resolution up to seconds

    ex:
        time2mag_grid_file([2003,11,20])
            returns l
        time2mag_grid_file([2003,11,20,7,0,11])
            returns 'mag_grid_e20031120-070011.out'
        time2mag_grid_file([2003,11,20,7,0,11,123])
            returns 'mag_grid_e20031120-070011.out'
    """
    time = list(time)
    filename = 'mag_grid_e' \
        + '%04d%02d%02d-%02d%02d%02d' % tpad(time, length=6) + '.out'
    return filename

def mag_grid_file2time(filename):
    """Extract time stamp from file name"""

    tstr = filename[10:] 
    y, m, d = int(tstr[0:4]), int(tstr[4:6]), int(tstr[6:8])
    h, M, s = int(tstr[9:11]), int(tstr[11:13]), int(tstr[13:15])
    f = 0
    return [y, m, d, h, M, s, f]

def getdata(run, time, debug=True):
    """

    Parameters
    ----------!!!!!

    Returns
    -------
    list of np arrays = [data, headers]
    
    where data is (N,17) array of N datapoints given in the file and 
    headers is (17,) array of string for the corresponding quantities of data
    
    if the file doesnt exist, it is downloaded form mag.gmu.edu
    
    """

    filename = conf[run + '_cdf'] + time2mag_grid_file(time)

    if not os.path.exists(filename):
        dlfile(filename)
        if debug:
            print('Downloading' + filename)

    data = np.genfromtxt(filename, skip_header=4)

    f = open(filename, 'r')
    first = f.readline()
    second = f.readline()
    third = f.readline()
    headerline = f.readline()
    f.close()

    headers = headerline[:-1].split(' ')
    assert(first[:19] == 'Magnetometer grid (')
    csyst = first[19:22]

    #headers = np.loadtxt(filename, dtype=str, skiprows=3, max_rows=1)

    return [data, headers, csyst]

def analyzedata(filename, MLAT, MLON, debug=True):
    """

    Parameters
    ----------
    filename : string or tuple/list/array
        if tuple/list/array, it's treated as time and converted to filename str
    MLAT : float
        input magnetic latitude in degrees.
    MLON : float
        input magnetic longitude in degrees.


    debug : boolean, OPTIONAL
    the default is True. 

    Returns
    -------
    ret : list 
        ret[0] is the index (int) of datapoint which the mlat and mlon value are near what was inputed
        ret[1],ret[2],ret[3] are the dBnMhd,dBeMhd,dBdMhd value at that datapoint
        
        if debug=True then this prints all the magnetic field information for that datapoint
        
    """
    data, headers = getdata(filename, debug=debug)
    if debug:
        print(headers)
    Tr = np.all([MLON-0.5 <= data[:, 0], 
                 data[:, 0] <= MLON+0.5, MLAT-0.5 <= data[:, 1],
                 data[:, 1] <= MLAT+0.5], axis=0)
    k = np.where(Tr==True)[0][0]
    if debug:
        print(k)
        print(data[k, 0],data[k, 1])
    ret = [k]
    for i in range(3): # i=0 <-> north , i=1 <-> east , i=3 <-> down
        if debug:
            print(headers[2+i], headers[5+i], headers[8+i], headers[11+i], headers[14+i], 'sum should equal full dB')
            print(data[k, 2+i], data[k, 5+i], data[k, 8+i], data[k, 11+i], data[k, 14+i], data[k, 5+i] + data[k, 8+i] + data[k, 11+i] + data[k, 14+i])
        ret.append(data[k, 5+i])
    dB_norm_Mhd = np.sqrt(data[k, 5+0]**2 + data[k, 5+1]**2 + data[k, 5+2]**2)
    if debug:
        print(headers[5+0] + '  ' + headers[5+1] + '  '+  headers[5+2])
        print('dB_norm_Mhd = ' + str(dB_norm_Mhd))

    return ret
