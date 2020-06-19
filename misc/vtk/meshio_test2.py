# meshio_test2

import sys
import numpy as np
import meshio

if sys.version_info.major != 2:
    raise Exception("Python 2 is required.")

in_fname = 'mesh_test2_import.vtk'
#out_fname = '/home/gary/magnetosphere/test_data/mesh_export.vtk'
cust_fname = 'meshio_test2_cust.vtk'

"""
Demonstrate writing VTK using meshio using data created in Python.
"""

a = 2
b = 2
c = 3

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
#B1, B2, B3 = np.meshgrid(X, Y, Z) # seems more natural but doesnt work with vtk structured_grid format
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

def unflat(k, npar=False):
    if k is None:
        return None
    if npar:
        return cell_inds[k,:]
    return tuple(cell_inds[k,:])

def flat(i,j,k):
    #https://stackoverflow.com/questions/16094563/numpy-get-index-where-value-is-true
    return np.where(np.all([cell_inds[:,0]==i,cell_inds[:,1]==j,cell_inds[:,2]==k],axis=0))[0][0]

def unflat2(k, npar=False):
    if k is None:
        return None
    if npar:
        return ind_p_np[k,:]
    return tuple(ind_p_np[k,:])

def flat2(i,j,k):
    #https://stackoverflow.com/questions/16094563/numpy-get-index-where-value-is-true
    return np.where(np.all([ind_p_np[:,0]==i,ind_p_np[:,1]==j,ind_p_np[:,2]==k],axis=0))[0][0]

if False:
    tup=(1,4,7)
    print(tup)
    print(flat(*tup))
    print(unflat(flat(*tup)))
    
    ind=15
    print(ind)
    print(unflat(ind))
    print(flat(*unflat(ind)))
    print(unflat(ind,npar=True))
    print(flat(*tuple(unflat(ind,npar=True))))


shifts = np.array( [[0, 0, 0],
                    [1, 0, 0],
                    [1, 1, 0],
                    [0, 1, 0],
                    [0, 0, 1],
                    [1, 0, 1],
                    [1, 1, 1],
                    [0, 1, 1]] )

print(shifts)

cust_cells = np.nan*np.empty(((a-1)*(b-1)*(c-1), 8))
for k in range((a-1)*(b-1)*(c-1)):
    for n in range(8):
        tup = tuple(unflat(k, npar=True) + shifts[n,:])
        cust_cells[k,n] = flat2(*tup)

# the structured grid is just a rectangular grid of points with rectangular prism cells in usual way.
# if there are a points along edge in x direction, then it will be made of a-1 segments of cells edges
# thus number of points is a*b*c and number of cells is (a-1)*(b-1)*(c-1)
print('points in x array = a =', a)
print('points in x array = b =', b)
print('points in x array = c =', c)
print('num points = a*b*c =', a*b*c)
print('num cells = (a-1)*(b-1)*(c-1) =', (a-1)*(b-1)*(c-1))



fil = open(in_fname,'w')

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

fil.close()

'''
fil = open(in_fname,'w')

fil.write('# vtk DataFile Version 2.0\n')
fil.write('Really cool data\n')
fil.write('BINARY\n')
fil.write('DATASET STRUCTURED_GRID\n')
fil.write('DIMENSIONS ' + str(a) + ' ' + str(b) + ' ' + str(c) + '\n')
fil.write('POINTS ' + str(a*b*c) + ' float\n')
fil.write(str(p_np[l,0]) + ' ' + str(p_np[l,1]) + ' ' + str(p_np[l,2]) + '\n')
fil.write('POINT_DATA ' + str(a*b*c) + '\n')
fil.write('SCALARS sample_scalars float 1\n')
fil.write('LOOKUP_TABLE my_table\n')
for l in range(p_np.shape[0]):
    fil.write(str(fvals[l]) + ' ')

fil.close()
'''

mesh = meshio.vtk.read(in_fname)

file_cells=mesh.cells_dict['hexahedron']
file_points=mesh.points

print(file_points.shape)
print(file_cells.shape)

print(np.all(p_np==file_points))
print(np.all(cust_cells==file_cells))

#print(cust_cells)
#print(cust_cells.astype(int))

c_np = [('hexahedron', cust_cells.astype(int))]
fvals_dict= {'sample_scalars': fvals}
cust_mesh = meshio.Mesh(p_np, c_np, point_data=fvals_dict)
meshio.vtk.write(cust_fname, cust_mesh, binary=False)
