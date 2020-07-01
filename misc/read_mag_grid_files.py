"""
This script reads the calculated values for the different contributions to the
magnetic field that are in mag_grid____.out files
"""

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
from urlretrieve import urlretrieve


def time2mag_grid_file(time):
    time = list(time)
    filename = 'mag_grid_e' \
        + '%04d%02d%02d-%02d%02d%02d' % tuple(time[0:6]) + '.out'
    return filename


def getdata(filename, MLON, MLAT, debug=True):
    if type(filename) != str:
        filename = time2mag_grid_file(filename)
    if not os.path.exists(conf['run_path'] + filename):
        urlretrieve(conf['run_url'] + filename, conf['run_path'] + filename)
        if debug:
            print('Downloading' + filename)
    fname = conf['run_path'] + filename

    data = np.genfromtxt(fname, skip_header=4)
    headers = np.loadtxt(fname, dtype=str, skiprows=3, max_rows=1)

    print(headers)

    Tr = np.all([MLON-0.5 <= data[:, 0], 
                 data[:, 0] <= MLON+0.5, MLAT-0.5 <= data[:, 1],
                 data[:, 1] <= MLAT+0.5], axis=0) # want 176.00, 57.50
    k = np.where(Tr==True)[0][0]
    print(k)
    print(data[k, 0],data[k, 1])
    i = 1
    print(headers[2+i], headers[5+i], headers[8+i], headers[11+i], headers[14+i], 'sum should equal full dB')
    print(data[k, 2+i], data[k, 5+i], data[k, 8+i], data[k, 11+i], data[k, 14+i], data[k, 5+i] + data[k, 8+i] + data[k, 11+i] + data[k, 14+i])