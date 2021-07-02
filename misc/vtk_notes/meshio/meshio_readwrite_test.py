"""
pip install meshio
only runs with python3

Typical output:
    
num points is  576000
1.70s: Write ASCII VTK with structured_grid using loop
1.23s: Read ASCII VTK with structured grid using meshio.read()
1.23s: Read ASCII VTK with structured grid using meshio.vtk.read()
1.22s: Read ASCII VTK with structured grid using meshio.vtk.read()
1.23s: Read ASCII VTK with structured grid using meshio.Mesh.read()
0.05s: Write with meshio.write()
0.04s: Write with meshio.vtk.write()
2.19s: Write with meshio.write(..., binary=False)
0.05s: Write with m.write(..., file_format="vtk")

Conclusions:
    1. Writing ASCII VTK is slightly faster using custom code (1.70 vs. 2.19)
    2. Writing binary using meshio is ~2.19/.05 ~ 50 times faster than writing
       ASCII with meshio.
    3. Best to use meshio.vtk.write(...,binary=True) since it can be more easily
       debuged by changing binary to False
"""

import os
import sys
import time
import numpy as np
import meshio

if sys.version_info.major != 3:
    raise Exception("Python 3 is required.")
      
in_fname = 'mesh_test1_import.vtk'
out_fname = 'mesh_test1_export.vtk'
if os.path.exists(in_fname): os.system(f'rm {in_fname}')

a = 60
b = 80
c = 120
print('Nx*Ny*Nz =', a*b*c)

X = np.array([i for i in range(a)])
Y = np.array([j for j in range(b)])
Z = np.array([k for k in range(c)])
B2, B3, B1 = np.meshgrid(Y, Z, X)
#B1, B2, B3 = np.meshgrid(X, Y, Z) # seems more natural but doesnt work with vtk structured_grid format
B1 = B1.flatten(order='C')
B2 = B2.flatten(order='C')
B3 = B3.flatten(order='C')
p_np = np.column_stack((B1, B2, B3))

fvals = 0.1*p_np[:,0]*p_np[:,1]*p_np[:,2]

to = time.time()
with open(in_fname,'w') as fil:
    fil.write('# vtk DataFile Version 2.0\n')
    fil.write('test data structured grid with scalar associated with each point\n')
    fil.write('ASCII\n')
    fil.write('DATASET STRUCTURED_GRID\n')
    fil.write('DIMENSIONS ' + str(a) + ' ' + str(b) + ' ' + str(c) + '\n')
    fil.write('POINTS ' + str(a*b*c) + ' int\n')
    for l in range(p_np.shape[0]):
        fil.write(str(p_np[l,0]) + ' ' + str(p_np[l,1]) + ' ' + str(p_np[l,2]) + '\n')
    fil.write('POINT_DATA ' + str(a*b*c) + '\n')
    fil.write('SCALARS sample_scalars float 1\n')
    fil.write('LOOKUP_TABLE my_table\n')
    for l in range(p_np.shape[0]):
        fil.write(str(fvals[l]) + ' ')

tf = time.time()
print('{0:.2f}s: Write ASCII VTK with structured_grid using loop'.format(tf - to))

def read_method(method):
    if method==0:
        to = time.time()
        mesh = meshio.read(
            in_fname,  # string, os.PathLike, or a buffer/open file
            file_format="vtk"  # optional if filename is a path; inferred from extension
        )
        # when using python 2, gives error: Unknown file format 'vtk' 
        # if the optional file_format isn't given, gives error: Only VTK
        # UNSTRUCTURED_GRID supported (not STRUCTURED_GRID)
        tf = time.time()
        print('{0:.2f}s: Read ASCII VTK with structured grid using meshio.read()'.format(tf - to))

    if method==1:
        to = time.time()
        mesh = meshio.vtk.read(in_fname)
        tf = time.time()
        print('{0:.2f}s: Read ASCII VTK with structured grid using meshio.vtk.read()'.format(tf - to))

    if method==2:
        to = time.time()
        mesh = meshio.Mesh.read(in_fname, "vtk")  # same arguments as meshio.read
        tf = time.time()
        print('{0:.2f}s: Read ASCII VTK with structured grid using meshio.Mesh.read()'.format(tf - to))

def write_method(method):
    mesh = meshio.vtk.read(in_fname)
    # all are comparable speed except meshio.vtk.write(...,binary=False)
    if method==0:
        to = time.time()
        meshio.write(
            out_fname,  # str, os.PathLike, or buffer/ open file
            mesh,
            file_format="vtk",  # optional if first argument is a path; inferred from extension
        )
        tf = time.time()
        print('{0:.2f}s: Write with meshio.write()'.format(tf - to))
    
    if method==1:
        to = time.time()
        # writes ascii unstructured grid but with cells
        meshio.vtk.write(out_fname, mesh, binary=True) 
        tf = time.time()
        print('{0:.2f}s: Write with meshio.vtk.write()'.format(tf - to))

    if method==2:
        to = time.time()
        meshio.vtk.write(out_fname, mesh, binary=False)
        tf = time.time()
        print('{0:.2f}s: Write with meshio.write(..., binary=False)'.format(tf - to))
    
    if method==3:
        to = time.time()
        # same arguments as meshio.write, besides `mesh`
        mesh.write(out_fname, file_format="vtk")  
        tf = time.time()
        print('{0:.2f}s: Write with mesh.write(..., file_format="vtk")'.format(tf - to))

for i in range(4):
    read_method(i)

for i in range(4):
    write_method(i)
    os.system(f'rm {out_fname}')
