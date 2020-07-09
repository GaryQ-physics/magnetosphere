import sys
import os
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../../' )
from config import conf

import spacepy.pybats.bats as bats

# read in the 3d magnetosphere
filename = conf['run_path'] + "3d__var_3_e20031120-070000-000.out"
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

Test = False

D = 3.96875

xax = np.linspace(-D, D, 128)
yax = np.linspace(-D, D, 128)
zax = np.linspace(-D, D, 128)

X = np.empty((0, 3)) # Combined slices in x
X2 = np.empty((0, 3)) # Combined slices in z
Tr = np.empty((0, ), dtype=bool)
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

#X2 = np.column_stack([x[Tr], y[Tr], z[Tr]])
full_grid = True

if Test:
    import biot_savart as bs
    X = bs.make_grid(xax, yax, zax, 0, 0, 0)[0]
if full_grid:
    X = np.column_stack([x, y, z])

print(X)
print(X.shape)
print(X2)
print(X2.shape)


list1 = list(X)
list2 = list(X2)

ltup1 = [tuple(l) for l in list1]
ltup2 = [tuple(l) for l in list2]

consistent = set(ltup1) == set(ltup2)

if consistent:
    print('consistent')
else:
    print('inconsistent')
if not (Test or full_grid):
    assert(consistent)

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
