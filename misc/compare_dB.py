"""
This script reads the calculated values for the different contributions to the magnetic field that are in mag_grid____.out files

The goal is to compare them with values calculated by biot-savart integral in this repository.
"""

link = 'http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/mag_grid_e20031120-070000.out'
#conf['run_url']=='http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/'

import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

debug = True

#import cut_plane_plot as cp
from urlretrieve import urlretrieve

if not os.path.exists(conf['run_path'] + 'mag_grid_e20031120-070000.out'):
    urlretrieve(conf['run_url'] + 'mag_grid_e20031120-070000.out', conf['run_path'] + 'mag_grid_e20031120-070000.out')

fname = conf['run_path'] + 'mag_grid_e20031120-070000.out'
#with open(conf['run_path'] + 'mag_grid_e20031120-070000.out','r') as f:
#    lines = f.readlines()

data = np.genfromtxt(fname, skip_header=4)
headers=np.loadtxt(fname, dtype=str, skiprows=3, max_rows=1)

print(headers)

Tr = np.all([176.00 == data[:, 0], 57. < data[:, 1], data[:, 1] < 58.], axis=0) # want 176.00, 57.50
k = np.where(Tr==True)[0][0]
print(k)
print(data[k, 0],data[k, 1])
i = 1
print(headers[2+i], headers[5+i], headers[8+i], headers[11+i], headers[14+i], 'sum should equal full dB')
print(data[k, 2+i], data[k, 5+i], data[k, 8+i], data[k, 11+i], data[k, 14+i], data[k, 5+i] + data[k, 8+i] + data[k, 11+i] + data[k, 14+i])
