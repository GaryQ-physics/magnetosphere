"""
pip install meshio

only runs with python3
"""

"""
Demonstrate creating a mesh object from numpy arrays created in Python
    and comparing to mesh object imported from existing vtk. The created
    mesh object can then be written to a new vtk. 

Everything is written explicitly without reference to external function, Including the making of the grids

Conclusion:
       writing the numpy data for a structured grid is possible using meshio, but the output file structure is not 
        an actual structured grid, instead it's an unstructured grid with cells. 
       Was not able to determine how to write structured grid using meshio.

In paraview 5.8.0, for the unstructured grid with cells file format, the "Slice With Plane" filter must be used instead of just "Slice"

"""

import sys
import numpy as np
import meshio

if sys.version_info.major != 3:
    raise Exception("Python 3 is required.")

fname1 = 'mesh_test2_ToImport.vtk'
fname2 = 'meshio_test2_FromNumpyArrays.vtk'

a = 20
b = 20
c = 30

############### Make Needed Grids
Xind = np.array([i for i in range(a)])
Yind = np.array([j for j in range(b)])
Zind = np.array([k for k in range(c)])
C2, C3, C1 = np.meshgrid(Yind, Zind, Xind)
C1 = C1.flatten(order='C')
C2 = C2.flatten(order='C')
C3 = C3.flatten(order='C')
ind_p_np = np.column_stack((C1, C2, C3))

X = np.array([1.5*i for i in range(a)])
Y = np.array([1.5*j for j in range(b)])
Z = np.array([1.5*k for k in range(c)])
B2, B3, B1 = np.meshgrid(Y, Z, X)
B1 = B1.flatten(order='C')
B2 = B2.flatten(order='C')
B3 = B3.flatten(order='C')
p_np = np.column_stack((B1, B2, B3))

fvals = 0.1*p_np[:,0]*p_np[:,1]*p_np[:,2]

Xc = np.array([i for i in range(a-1)])
Yc = np.array([j for j in range(b-1)])
Zc = np.array([k for k in range(c-1)])
D2, D3, D1 = np.meshgrid(Yc, Zc, Xc)
D1 = D1.flatten(order='C')
D2 = D2.flatten(order='C')
D3 = D3.flatten(order='C')
cell_inds = np.column_stack((D1, D2, D3))

# the structured grid is just a rectangular grid of points with rectangular prism cells in usual way.
# if there are a points along edge in x direction, then it will be made of a-1 segments of cells edges
# thus number of points is a*b*c and number of cells is (a-1)*(b-1)*(c-1)
print('points in x array = a =', a)
print('points in x array = b =', b)
print('points in x array = c =', c)
print('num points = a*b*c =', a*b*c)
print('num cells = (a-1)*(b-1)*(c-1) =', (a-1)*(b-1)*(c-1))


############### write ascii structured grid VTK grid via python loop
with open(fname1,'w') as fil:
    fil.write('# vtk DataFile Version 2.0\n')
    fil.write('Really cool data\n')
    fil.write('ASCII\n')
    fil.write('DATASET STRUCTURED_GRID\n')
    fil.write('DIMENSIONS ' + str(a) + ' ' + str(b) + ' ' + str(c) + '\n')
    fil.write('POINTS ' + str(a*b*c) + ' float\n')
    for l in range(p_np.shape[0]):
        fil.write(str(p_np[l,0]) + ' ' + str(p_np[l,1]) + ' ' + str(p_np[l,2]) + '\n')
    fil.write('POINT_DATA ' + str(a*b*c) + '\n')
    fil.write('SCALARS sample_scalars float 1\n')
    fil.write('LOOKUP_TABLE my_table\n')
    for l in range(p_np.shape[0]):
        fil.write(str(fvals[l]) + ' ')

############### Import mesh object (file_mesh) from the above written vtk
file_mesh = meshio.vtk.read(fname1)

file_cells=file_mesh.cells_dict['hexahedron']
file_points=file_mesh.points

print(file_points.shape)
print(file_cells.shape)


############### Create custom mesh object (cust_mesh) from numpy arrays
shifts = np.array( [[0, 0, 0],
                    [1, 0, 0],
                    [1, 1, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                    [1, 0, 1],
                    [1, 1, 1],
                    [0, 1, 1]] )

cust_cells = np.nan*np.empty(((a-1)*(b-1)*(c-1), 8))
for k in range((a-1)*(b-1)*(c-1)):
    for n in range(8):
        tup = tuple(cell_inds[k,:] + shifts[n,:])

        cust_cells[k,n] = np.where(np.all([ind_p_np[:,0]==tup[0],ind_p_np[:,1]==tup[1],ind_p_np[:,2]==tup[2]],axis=0))[0][0]

c_np = [('hexahedron', cust_cells.astype(int))]
fvals_dict= {'sample_scalars': fvals}
cust_mesh = meshio.Mesh(p_np, c_np, point_data=fvals_dict)
meshio.vtk.write(fname2, cust_mesh, binary=False)

# Compair to make sure atributes are the same
##################
print(np.all(file_points==p_np))
print(np.all(file_mesh.points==cust_mesh.points))
print(np.all(file_cells==cust_cells))
print(np.all(file_mesh.cells_dict['hexahedron']==cust_mesh.cells_dict['hexahedron']))
