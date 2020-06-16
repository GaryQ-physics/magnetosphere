"""
pip install meshio

only runs with python3
"""

# meshio_test1

import numpy as np
import meshio
#import pathlib
import time

# Older Version?  https://python.hotexamples.com/examples/meshio/-/write/python-write-function-examples.html


in_fname = '/home/gary/magnetosphere/test_data/mesh_import.vtk'
out_fname = '/home/gary/magnetosphere/test_data/mesh_export.vtk'


a = 60
b = 80
c = 120

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

print('num points is ', a*b*c)


to = time.time()


fil = open(in_fname,'w')

fil.write('# vtk DataFile Version 2.0\n')
fil.write('Really cool data\n')
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

fil.close()

tf = time.time()
print('time to write structured_grid with loop =', tf - to)

#print(mesh1.points)
#print(mesh1.cells)

Try=2

if Try==1:
    to = time.time()
    mesh1 = meshio.read(
        in_fname,  # string, os.PathLike, or a buffer/open file
        file_format="vtk"  # optional if filename is a path; inferred from extension
    )
    """
    when using python 2, gives error: Unknown file format 'vtk' 
        if the optional file_format isn't given, gives error: Only VTK UNSTRUCTURED_GRID supported (not STRUCTURED_GRID)
    """

    tf = time.time()
    print('time to read the ascii structured grid with meshio.read( ) =', tf - to)

if Try==2 or Try==3:
    to = time.time()
    mesh2 = meshio.vtk.read(in_fname)
    tf = time.time()
    print('time to read the ascii structured grid with meshio.vtk.read( ) =', tf - to)

if Try==4:
    to = time.time()
    m = meshio.Mesh.read(in_fname, "vtk")  # same arguments as meshio.read
    tf = time.time()
    print('time to read the ascii structured grid with meshio.Mesh.read( ) =', tf - to)

#print(mesh1==mesh2)


if Try==1:
    to = time.time()
    meshio.write(
        out_fname,  # str, os.PathLike, or buffer/ open file
        mesh1,
        file_format="vtk",  # optional if first argument is a path; inferred from extension
    )
    tf = time.time()
    print('time to write with meshio.write( ) =', tf - to)

if Try==2:
    to = time.time()
    meshio.vtk.write(out_fname, mesh2, binary=True) #writes ascii unstructured grid but with cells
    tf = time.time()
    print('time to write with meshio.vtk.write( ,binary=True) =', tf - to)

if Try==3:
    to = time.time()
    meshio.vtk.write(out_fname, mesh2, binary=False)
    tf = time.time()
    print('time to write with meshio.vtk.write( ,binary=False) =', tf - to)

if Try==4:
    to = time.time()
    m.write(out_fname, file_format="vtk")  # same arguments as meshio.write, besides `mesh`
    tf = time.time()
    print('time to write with m.write( ,binary) =', tf - to)

# all are comparable speed except meshio.vtk.write( ,binary=False)
# best to use meshio.vtk.write( ,binary=True) since it can be more easily debuged by changing binary to False

if False:
    points = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        ])
    cells = [
        ("triangle", np.array([[0, 1, 2]]))
    ]
    '''
    meshio.write_points_cells(
        out_fname,
        points,
        cells,
        # Optionally provide extra data on points, cells, etc.
        # point_data=point_data,
        # cell_data=cell_data,
        # field_data=field_data
        )
    '''
    
    mesh = meshio.Mesh(points, cells)
    meshio.write(
        out_fname,  # str, os.PathLike, or buffer/ open file
        mesh,
        file_format="vtk",  # optional if first argument is a path; inferred from extension
    )
    print(mesh.points)