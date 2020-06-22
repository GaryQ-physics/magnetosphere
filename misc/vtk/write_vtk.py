# write_vtk.py - Demonstrate several methods of writing a VTK file
# 1. ASCII
# 2. Binary using loop
# 3. Binary without loop

# To read files in Paraview (tested on 5.4.1 and 5.8.0):
# File -> Open one of the created VTK files
# Click Apply
# Under Display(Geometry Representation), select Glyph
# Click Apply

# Works in Python 2 and 3  (NOT FOR ME!!!!!!!!!!)

import numpy as np

N = 3
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
f = open(fname,'wb')
f.write('# vtk DataFile Version 3.0\n'.encode('ascii'))
f.write('Wind velocity at unstructured points\n'.encode('ascii'))
f.write('BINARY\n'.encode('ascii'))
f.write('\n'.encode('ascii'))
f.write('DATASET POLYDATA\n'.encode('ascii'))
s = 'POINTS '+str(N**3)+' float\n'
f.write(s.encode('ascii'))
for i in range(N):
    for j in range(N):
        for k in range(N):
            f.write(struct.pack('>f', float(X[i])))
            f.write(struct.pack('>f', float(Y[j])))
            f.write(struct.pack('>f', float(Z[k])))

f.write('\n'.encode('ascii'))
s = 'POINT_DATA '+str(N**3)+'\n'
f.write(s.encode('ascii'))
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
f = open(fname,'wb'.encode('ascii'))
f.write('# vtk DataFile Version 3.0\n'.encode('ascii'))
f.write('Wind velocity at unstructured points\n'.encode('ascii'))
f.write('BINARY\n'.encode('ascii'))
f.write('\n'.encode('ascii'))
f.write('DATASET POLYDATA\n'.encode('ascii'))
s = 'POINTS '+str(N**3)+' float\n'
f.write(s.encode('ascii'))
f.write(P.tobytes())

f.write('\n'.encode('ascii'))
s = 'POINT_DATA '+str(N**3)+'\n'
f.write(s)
f.write('VECTORS point_vectors float\n'.encode('ascii'))
f.write(P.tobytes())

f.close()

###############################################################################
