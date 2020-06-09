link = 'http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/mag_grid_e20031120-070000.out'
#conf['run_url']=='http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/'

import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

debug = True

import cut_plane_plot as cp
from urlretrieve import urlretrieve

if not os.path.exists(conf['run_path'] + 'mag_grid_e20031120-070000.out'):
    urlretrieve(conf['run_url'] + 'mag_grid_e20031120-070000.out', conf['run_path'] + 'mag_grid_e20031120-070000.out')

fname = conf['run_path'] + 'mag_grid_e20031120-070000.out'
#with open(conf['run_path'] + 'mag_grid_e20031120-070000.out','r') as f:
#    lines = f.readlines()

'''
f = open(conf['run_path'] + 'mag_grid_e20031120-070000.out', 'r')
lines = f.readlines()

headers = [str(j) for j in lines[3].split()] 
print(headers)

N = len(lines)-4
#N=10
for i in range(N):
    k = i+4
    data = [float(j) for j in lines[k].split()] #numpy.genfromtxt .  what about numpy.loadtxt?
    if 176.00 == data[0] and 57. < data[1] and data[1] < 58.: # want 176.00, 57.50
        print(k)
        print(data[0],data[1])
        print(headers[2], headers[5], headers[8], headers[11], headers[14])
        print(data[2], data[5], data[8], data[11], data[14], data[5] + data[8] + data[11] + data[14])
'''

data = np.genfromtxt(fname, skip_header=4)
headers=np.loadtxt(fname, dtype=str, skiprows=3, max_rows=1)
'''
print(headers)
for k in range(data.shape[0]):
    if 176.00 == data[k, 0] and 57. < data[k, 1] and data[k, 1] < 58.: # want 176.00, 57.50
        print(k)
        print(data[k, 0],data[k, 1])
        print(headers[2], headers[5], headers[8], headers[11], headers[14])
        print(data[k, 2], data[k, 5], data[k, 8], data[k, 11], data[k, 14], data[k, 5] + data[k, 8] + data[k, 11] + data[k, 14])
'''

Tr = np.all([176.00 == data[:, 0], 57. < data[:, 1], data[:, 1] < 58.], axis=0)# want 176.00, 57.50
k = np.where(Tr==True)[0][0]
print(k)
print(data[k, 0],data[k, 1])
print(headers[2], headers[5], headers[8], headers[11], headers[14])
print(data[k, 2], data[k, 5], data[k, 8], data[k, 11], data[k, 14], data[k, 5] + data[k, 8] + data[k, 11] + data[k, 14])
