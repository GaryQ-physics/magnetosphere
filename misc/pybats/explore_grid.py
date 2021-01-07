import sys
import os
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../../' )

from config import conf
import util

import spacepy.pybats.bats as bats

# read in the 3d magnetosphere
#filename = conf['run_path'] + "3d__var_3_e20031120-070000-000.out"
#filename = conf['SWPC_raw'] + "3d__var_1_t00000000_n0002500.out"
filename = util.time2CDFfilename('DIPTSUR2',(2019,9,2,6,30,0))[:-4]
#filename = conf['SCARR5'+'_cdf'] + "3d__var_3_e20031120-070000-000.out"
data3d = bats.Bats2d(filename)



import numpy as np
# get the cell coordinates
x = data3d['x']
y = data3d['y']
z = data3d['z']
#bz = data3d['bz']
#jy = data3d['jy']
#jz = data3d['jz']
#p = data3d['p']

x = np.array(x, dtype=float)
y = np.array(y, dtype=float)
z = np.array(z, dtype=float)
#p = np.array(p)


D = 3.96875

xax = np.linspace(-D, D, 128)
yax = np.linspace(-D, D, 128)
zax = np.linspace(-D, D, 128)

X = np.empty((0, 3)) # Combined slices in x
X2 = np.empty((0, 3)) # Combined slices in z
for i in range(128):
    tr = x == i*0.0625 - D

    x_ = x[tr]
    y_ = y[tr]
    z_ = z[tr]

    S = np.column_stack([x_, y_, z_]) # each slice
    X = np.vstack((X, S))
    #Tr = np.hstack((Tr,tr))

    tr2 = z == i*0.0625 - D

    x_2 = x[tr2]
    y_2 = y[tr2]
    z_2 = z[tr2]

    S2 = np.column_stack([x_2, y_2, z_2]) # each slice
    X2 = np.vstack((X2, S2))

Tr = np.all([-D<=x, x<=D, -D<=y, y<=D, -D<=z, z<=D], axis=0)
x_3 = x[Tr]
y_3 = y[Tr]
z_3 = z[Tr]
X3 = np.column_stack([x_3, y_3, z_3])


Test = False
full_grid = False
if Test:
    import biot_savart as bs
    X = bs.make_grid(xax, yax, zax, 0, 0, 0)[0]
if full_grid:
    X = np.column_stack([x, y, z])


print(X)
print(X.shape)
print(X2)
print(X2.shape)
print(X3)
print(X3.shape)

list1 = list(X)
list2 = list(X2)
list3 = list(X3)

ltup1 = [tuple(l) for l in list1]
ltup2 = [tuple(l) for l in list2]
ltup3 = [tuple(l) for l in list3]

set1 = set(ltup1)
set2 = set(ltup2)
set3 = set(ltup3)

consistent = set1 == set2

if consistent:
    print('consistent')
else:
    print('inconsistent')

print('check edges: ' + str(set1 == set3))

#if not (Test or full_grid):
#    assert(consistent)

difs = np.array(list(set3.difference(set1)))
Rdifs = np.sqrt(difs[:,0]**2 + difs[:,1]**2 + difs[:,2]**2)
print(np.min(Rdifs))
print(Rdifs.shape)

if Test:
    fname = 'explore_grid_linspace.vtk'
else:
    fname = 'explore_grid_sp.vtk'
if full_grid:
    fname = 'explore_grid_full.vtk'

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


if False:
    #print(x_)
    #print(y_)
    #print(z_)
    print(np.min(z_), np.max(z_))
    print(z_.shape)

    #print(yax)

    Y,Z = np.meshgrid(yax, zax)

    import matplotlib.pyplot as plt
    plt.plot(y_, z_, marker='.', color='b', linestyle='')
    plt.plot(Y.flatten(), Z.flatten(), marker='.', color='r', linestyle='', markersize=2.)
    plt.xlabel('xlable')
    plt.ylabel('ylable')
    #plt.axvline(x=1.25)
    #plt.axvline(x=1.2)
    plt.show()
