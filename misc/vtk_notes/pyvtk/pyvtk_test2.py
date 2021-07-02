"""
pip install pyvtk
runs with python2 or python3

Compare writing basic structured grid using pyvtk and looping f.write

adjust loop to determing wether you use loop or pyvtk

Typical output (python2):
    loop = True  time = 5.6
    loop = False  time = 1.5
Typical output (python3):
    loop = True  time = 3.6
    loop = False  time = 1.8

Conclusion:
    looping and writing the ascii is 4 (2) times faster than pyvtk even running binary in python 2 (3) , so pyvtk is prohibitively slow
"""

#import sys
#sys.path = ['..']+sys.path
import numpy as np
import pyvtk
import time

loop = True


a = 60
b = 80
c = 120

def f(x,y,z):
    return x*y*z


# -- ----------
X = np.array([i for i in range(a)])
Y = np.array([j for j in range(b)])
Z = np.array([k for k in range(c)])
B2, B3, B1 = np.meshgrid(Y, Z, X)
#B1, B2, B3 = np.meshgrid(X, Y, Z) # seems more natural but doesnt work with vtk structured_grid format
B1 = B1.flatten(order='C')
B2 = B2.flatten(order='C')
B3 = B3.flatten(order='C')
p_np = np.column_stack((B1, B2, B3))

fvals = f(p_np[:,0], p_np[:,1], p_np[:,2])

fname = 'pyvtk_test3'

#------------------

to = time.time()


if loop:
    fil = open(fname + '.vtk','w')

    fil.write('# vtk DataFile Version 2.0\n')
    fil.write('Really cool data\n')
    fil.write('ASCII\n')
    fil.write('DATASET STRUCTURED_GRID\n')
    fil.write('DIMENSIONS ' + str(a) + ' ' + str(b) + ' ' + str(c) + '\n')
    fil.write('POINTS ' + str(a*b*c) + ' int\n')
    for l in range(p_np.shape[0]):
        fil.write(str(p_np[l,0]) + ' ' + str(p_np[l,1]) + ' ' + str(p_np[l,2]) + '\n')
    fil.write('POINT_DATA ' + str(a*b*c) + '\n')
    fil.write('SCALARS sample_scalars int 1\n')
    fil.write('LOOKUP_TABLE my_table\n')
    for l in range(p_np.shape[0]):
        fil.write(str(fvals[l]) + ' ')
    fil.close()
else:
    pointdata = pyvtk.PointData(\
    pyvtk.Scalars(fvals,
            name='sample_scalars',
            lookup_table='my_table'),
            )
    vtk = pyvtk.VtkData(pyvtk.StructuredGrid([a,b,c],p_np), pointdata)
    vtk.tofile(fname,'binary')
    #vtk.tofile(fname,'ascii')


tf = time.time()

print('--------------')
print('num points is ', a*b*c)
print('time = ', tf-to)
print('loop = ', loop)
print('--------------')
