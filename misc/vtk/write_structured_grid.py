# write_structured_grid.py - Demonstrate methods of writing a structured grid VTK file
# 1. ASCII
# 3. Binary

# To read files in Paraview (tested on 5.8.0):
# File -> Open one of the created VTK files
# Click Apply
# select slice button
# Click Apply

# works in python2, but not python3 where it gives error:
#       write() argument must be str, not bytes

import numpy as np

var = 'Test'
N = 100
X=np.linspace(-5.,5.,N)
Y=np.linspace(-5.,5.,N)
Z=np.linspace(-5.,5.,N)

B2, B3, B1 = np.meshgrid(Y, Z, X)
B1 = B1.flatten(order='C')
B2 = B2.flatten(order='C')
B3 = B3.flatten(order='C')
Bgrid = np.column_stack((B1, B2, B3))

data = Bgrid[:,0]**2 + Bgrid[:,1]**2 + Bgrid[:,2]**2

###############################################################################

fname = 'structured_grid_test_ascii.vtk';

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Structured Grid ' + var + '\n')
f.write('ASCII\n')
f.write('DATASET STRUCTURED_GRID\n')
f.write('DIMENSIONS ' + str(N) + ' ' + str(N) + ' ' + str(N) + '\n' )
f.write('POINTS '+str(N*N*N)+' float\n')

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
f.write('POINT_DATA ' + str(N*N*N) + '\n')
f.write('SCALARS ' + var + ' float 1\n')
f.write('LOOKUP_TABLE default\n')
for l in range(Bgrid.shape[0]):
    f.write('%e\n' % (data[l],))

f.close()
print("Wrote " + fname)

###############################################################################

fname = 'structured_grid_test_binary.vtk';

Bgrid_toB = np.array(Bgrid, dtype='>f')
data_toB = np.array(data, dtype='>f')

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Structured Grid ' + var + '\n')
f.write('BINARY\n')
f.write('DATASET STRUCTURED_GRID\n')
f.write('DIMENSIONS ' + str(N) + ' ' + str(N) + ' ' + str(N) + '\n' )
f.write('POINTS '+str(N*N*N)+' float\n')

f.write(Bgrid_toB.tobytes())

f.write('\n')
f.write('POINT_DATA ' + str(N*N*N) + '\n')
f.write('SCALARS ' + var + ' float 1\n')
f.write('LOOKUP_TABLE default\n')

f.write(data_toB.tobytes())

f.close()
print("Wrote " + fname)

###############################################################################
