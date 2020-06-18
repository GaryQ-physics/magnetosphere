# write_vtk.py - Demonstrate several methods of writing a VTK file
# 1. ASCII
# 2. Binary using loop
# 3. Binary without loop

# two different conventions for enumerating the grid are shown

# To read files in Paraview (tested on 5.4.1 and 5.8.0):
# File -> Open one of the created VTK files
# Click Apply
# Under Display(Geometry Representation), select Glyph
# Click Apply

# works in python2, but not python3 where it gives error:
#       write() argument must be str, not bytes

import numpy as np

N = 100
X=np.linspace(-5.,5.,N)
Y=np.linspace(-5.,5.,N)
Z=np.linspace(-5.,5.,N)

###############################################################################

fname = 'write_vtk_ascii.vtk';

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Wind velocity at unstructured points\n')
f.write('ASCII\n')
f.write('\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS '+str(N**3)+' float\n')
for i in range(N):
    for j in range(N):
        for k in range(N):
            f.write('%e %e %e\n' % (X[i],Y[j],Z[k]))

f.write('\n')
f.write('POINT_DATA '+str(N**3)+'\n')
f.write('VECTORS point_vectors float\n')
for i in range(N):
    for j in range(N):
        for k in range(N):
            f.write('%e %e %e\n' % (X[i],Y[j],Z[k]))

f.close()

###############################################################################

import struct

fname = 'write_vtk_binary.vtk';

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Wind velocity at unstructured points\n')
f.write('BINARY\n')
f.write('\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS '+str(N**3)+' float\n')
for i in range(N):
    for j in range(N):
        for k in range(N):
            f.write(struct.pack('>f', float(X[i])))
            f.write(struct.pack('>f', float(Y[j])))
            f.write(struct.pack('>f', float(Z[k])))

f.write('\n')
f.write('POINT_DATA '+str(N**3)+'\n')
f.write('VECTORS point_vectors float\n')
for i in range(N):
    for j in range(N):
        for k in range(N):
            f.write(struct.pack('>f', float(X[i])))
            f.write(struct.pack('>f', float(Y[j])))
            f.write(struct.pack('>f', float(Z[k])))

f.close()

###############################################################################

fname = 'write_vtk_binary_noloop.vtk';

yv, xv, zv = np.meshgrid(X, Y, Z)
P = np.array([xv.flatten(), yv.flatten(), zv.flatten()], order='F')
P = P.flatten(order='F')
P = np.array(P, dtype='>f')

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Wind velocity at unstructured points\n')
f.write('BINARY\n')
f.write('\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS '+str(N**3)+' float\n')
f.write(P.tobytes())

f.write('\n')
f.write('POINT_DATA '+str(N**3)+'\n')
f.write('VECTORS point_vectors float\n')
f.write(P.tobytes())

f.close()

###############################################################################

'''
if;
yv, xv, zv = np.meshgrid(X, Y, Z)
P = np.array([xv.flatten(), yv.flatten(), zv.flatten()], order='F') , P.shape == (3,N)

B2, B3, B1 = np.meshgrid(Y, Z, X)
B1 = B1.flatten(order='C')
B2 = B2.flatten(order='C')
B3 = B3.flatten(order='C')
Bgrid = np.column_stack((B1, B2, B3)) , Bgrid.shape == (N,3)

then;
xv.flatten(order='F') == B1.flatten(order='C')
xv.flatten(order='C') == B1.flatten(order='F')
.flatten() defaults  order='C'
Bgrid.flattens and P.flattens arent equal in any combo
'''

fname = 'write_vtk_ascii_bsgrid.vtk';

B2, B3, B1 = np.meshgrid(Y, Z, X)
B1 = B1.flatten(order='C')
B2 = B2.flatten(order='C')
B3 = B3.flatten(order='C')
Bgrid = np.column_stack((B1, B2, B3))

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Wind velocity at unstructured points\n')
f.write('ASCII\n')
f.write('\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS '+str(N**3)+' float\n')

for l in range(Bgrid.shape[0]):
    f.write('%e %e %e\n' % (Bgrid[l,0], Bgrid[l,1], Bgrid[l,2]))

'''
should be equivalent to:
for k in range(N):
    for j in range(N):
        for i in range(N):
            f.write('%e %e %e\n' % (X[i],Y[j],Z[k]))

'''

f.write('\n')
f.write('POINT_DATA '+str(N**3)+'\n')
f.write('VECTORS point_vectors float\n')

for k in range(N):
    for j in range(N):
        for i in range(N):
            f.write('%e %e %e\n' % (X[i],Y[j],Z[k]))

'''
should be equivalent to:
for l in range(Bgrid.shape[0]):
    f.write('%e %e %e\n' % (Bgrid[l,0], Bgrid[l,1], Bgrid[l,2]))
'''

f.close()

###############################################################################

fname = 'write_vtk_binary_bsgrid.vtk';

B2, B3, B1 = np.meshgrid(Y, Z, X)
B1 = B1.flatten(order='C')
B2 = B2.flatten(order='C')
B3 = B3.flatten(order='C')
Bgrid = np.column_stack((B1, B2, B3))
Bgrid = np.array(Bgrid, dtype='>f')

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Wind velocity at unstructured points\n')
f.write('Binary\n')
f.write('\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS '+str(N**3)+' float\n')
f.write(Bgrid.tobytes())

f.write('\n')
f.write('POINT_DATA '+str(N**3)+'\n')
f.write('VECTORS point_vectors float\n')
f.write(Bgrid.tobytes())

f.close()

###############################################################################
