# pyvista only works with python3
import numpy as np

vectordata = False
N=10

xax = yax = zax = 1.*np.arange(N)

z = np.kron(np.ones(N),zax)
y = np.kron(yax,np.ones(N))
x = np.kron(xax, np.ones(N**2))
z = np.kron(np.ones(N), z)
y = np.kron(np.ones(N), y)

P = np.column_stack([x,y,z])

ind = np.arange(N**3).reshape(N,N,N)

E = []
for i in range(N-1):
    for j in range(N):
        for k in range(N):
            E.append((ind[i,j,k], ind[i+1,j,k]))
for i in range(N):
    for j in range(N-1):
        for k in range(N):
            E.append((ind[i,j,k], ind[i,j+1,k]))
for i in range(N):
    for j in range(N):
        for k in range(N-1):
            E.append((ind[i,j,k], ind[i,j,k+1]))
E = np.array(E, dtype=int)

F = []
for i in range(N-1):
    for j in range(N-1):
        for k in range(N):
            F.append( (ind[i,j,k], ind[i+1,j,k], ind[i+1,j+1,k], ind[i,j+1,k]) )
for i in range(N-1):
    for j in range(N):
        for k in range(N-1):
            F.append( (ind[i,j,k], ind[i+1,j,k], ind[i+1,j,k+1], ind[i,j,k+1]) )
for i in range(N):
    for j in range(N-1):
        for k in range(N-1):
            F.append( (ind[i,j,k], ind[i,j+1,k], ind[i,j+1,k+1], ind[i,j,k+1]) )
F = np.array(F, dtype=int)

V = []
for i in range(N-1):
    for j in range(N-1):
        for k in range(N-1):
            V.append( (ind[i,j,k], ind[i+1,j,k], ind[i+1,j+1,k], ind[i,j+1,k],
                       ind[i,j,k+1], ind[i+1,j,k+1], ind[i+1,j+1,k+1], ind[i,j+1,k+1])
                    )

V = np.array(V, dtype=int)



if False:
    perm = np.random.permutation(N**3)
    E = E[perm,:] #WRONG
    P = P[perm,:]

    print(E)
    print(F)
    print(P)

field = np.column_stack([ P[:,0]*P[:,2]**2+0.1 , P[:,0]*P[:,1]**2+0.1 , P[:,0]*P[:,2]+P[:,2]+0.1 ])
nE = E.shape[0]
nF = F.shape[0]
nV = V.shape[0]

'''
f = open('UNSTRUCTURED_GRID-volumesPeriodic.vtk','w')

f.write('# vtk DataFile Version 3.0\n')
f.write('Unstructured_grid cells\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('POINTS %d float\n'%(N**3))
np.savetxt(f,P, fmt = '%.3f')
f.write('CELLS %d %d\n'%(nV_Periodic, 9*nV_Periodic))
np.savetxt(f,np.column_stack([8*np.ones(nV_Periodic), V_Periodic]), fmt='%d')
f.write('CELL_TYPES %d\n'%(nV_Periodic))
np.savetxt(f, 12*np.ones(nV_Periodic), fmt = '%d') #12 here stands for hexahedron, see http://www.computingforscientists.info/ParaView

if vectordata:
    f.write('POINT_DATA %d\n'%(N**3))
    f.write('VECTORS field float\n')
    np.savetxt(f, field, fmt = '%.3f')
else:
    f.write('POINT_DATA %d\n'%(N**3))
    f.write('SCALARS magfield float 1\n')
    f.write('LOOKUP_TABLE default\n')
    np.savetxt(f, np.sqrt(np.einsum('ij,ij->i',field,field)), fmt = '%.3f')

f.close()
del f
'''


##### using V can generate streamlines ####
f = open('UNSTRUCTURED_GRID-volumes.vtk','w')

f.write('# vtk DataFile Version 3.0\n')
f.write('Unstructured_grid cells\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('POINTS %d float\n'%(N**3))
np.savetxt(f,P, fmt = '%.3f')
f.write('CELLS %d %d\n'%(nV,9*nV))
np.savetxt(f,np.column_stack([8*np.ones(nV), V]), fmt='%d')
f.write('CELL_TYPES %d\n'%(nV))
np.savetxt(f, 12*np.ones(nV), fmt = '%d') #12 here stands for hexahedron, see http://www.computingforscientists.info/ParaView

if vectordata:
    f.write('POINT_DATA %d\n'%(N**3))
    f.write('VECTORS field float\n')
    np.savetxt(f, field, fmt = '%.3f')
else:
    f.write('POINT_DATA %d\n'%(N**3))
    f.write('SCALARS magfield float 1\n')
    f.write('LOOKUP_TABLE default\n')
    np.savetxt(f, np.sqrt(np.einsum('ij,ij->i',field,field)), fmt = '%.3f')

f.close()
del f

if False:
    import pyvista as pv
    import vtk
    #pyvista outbuts "# vtk DataFile Version 5.1"  instead of 3.0 , and it has extra sections 
    cells = np.column_stack([8*np.ones(nV,dtype=int), V]).flatten()
    #offset = (8+1)*np.arange(nV)
    offset = np.zeros(1, dtype=int) # offset is a complete dummy here
    celltypes = np.empty(nV, dtype=np.uint8)
    celltypes[:] = vtk.VTK_HEXAHEDRON

    print('###dtypes:')
    print(offset.dtype)
    print(cells.dtype)
    print(celltypes.dtype)
    print(P.dtype)
    print('###\n\n')

    grid = pv.UnstructuredGrid(offset, cells, celltypes, P)

    print('\n\n')
    print('dummy inputed offset (offset):')
    print(offset)
    print('UnstructuredGrid computed offset (grid.offset):')
    print(grid.offset)

    grid.save('UNSTRUCTURED_GRID-volumes-pyvista.vtk', binary=False)


    ''' copy into python in terminal. will be the same for g and g2
    import pyvista as pv
    import numpy as np
    g = pv.read("UNSTRUCTURED_GRID-volumes-pyvista.vtk")
    g2 = pv.read("UNSTRUCTURED_GRID-volumes.vtk")
    type(g)
    #<class 'pyvista.core.pointset.UnstructuredGrid'>
    type(g2)
    #<class 'pyvista.core.pointset.UnstructuredGrid'>

    g.points
    g.cells
    g.offset

    g2.points
    g2.cells
    g2.offset
    ''' 



########################################################################
########################################################################

##### using E cannot generate streamlines ####
f = open('UNSTRUCTURED_GRID-edges.vtk','w')

f.write('# vtk DataFile Version 3.0\n')
f.write('Unstructured_grid cells\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('POINTS %d float\n'%(N**3))
np.savetxt(f,P, fmt = '%.3f')
f.write('CELLS %d %d\n'%(nE,3*nE))
np.savetxt(f,np.column_stack([2*np.ones(nE), E]), fmt='%d')
f.write('CELL_TYPES %d\n'%(nE))
np.savetxt(f, 3*np.ones(nE), fmt = '%d')

if vectordata:
    f.write('POINT_DATA %d\n'%(N**3))
    f.write('VECTORS field float\n')
    np.savetxt(f, field, fmt = '%.3f')
else:
    f.write('POINT_DATA %d\n'%(N**3))
    f.write('SCALARS magfield float 1\n')
    f.write('LOOKUP_TABLE default\n')
    np.savetxt(f, np.sqrt(np.einsum('ij,ij->i',field,field)), fmt = '%.3f')


f.close()
del f
##### using F cannot generate streamlines ####
f = open('UNSTRUCTURED_GRID-faces.vtk','w')

f.write('# vtk DataFile Version 3.0\n')
f.write('Unstructured_grid cells\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('POINTS %d float\n'%(N**3))
np.savetxt(f,P, fmt = '%.3f')
f.write('CELLS %d %d\n'%(nF,5*nF))
np.savetxt(f,np.column_stack([4*np.ones(nF), F]), fmt='%d')
f.write('CELL_TYPES %d\n'%(nF))
np.savetxt(f, 7*np.ones(nF), fmt = '%d')

if vectordata:
    f.write('POINT_DATA %d\n'%(N**3))
    f.write('VECTORS field float\n')
    np.savetxt(f, field, fmt = '%.3f')
else:
    f.write('POINT_DATA %d\n'%(N**3))
    f.write('SCALARS magfield float 1\n')
    f.write('LOOKUP_TABLE default\n')
    np.savetxt(f, np.sqrt(np.einsum('ij,ij->i',field,field)), fmt = '%.3f')


f.close()
del f

########################################################################
########################################################################



########################################################################
########################################################################
#cannot generate streamlines
f = open('POLYDATA-LINES.vtk','w')

f.write('# vtk DataFile Version 3.0\n')
f.write('Polydata Lines\n')
f.write('ASCII\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS %d float\n'%(N**3))
np.savetxt(f,P, fmt = '%.3f')
f.write('LINES %d %d\n'%(nE,3*nE))
np.savetxt(f, np.column_stack([2*np.ones(nE), E]), fmt='%d')

f.write('POINT_DATA %d\n'%(N**3))
#f.write('VECTORS field float\n')
np.savetxt(f, field, fmt = '%.3f')

f.close()
del f
#cannot generate streamlines
f = open('POLYDATA-POLYGONS.vtk','w')

f.write('# vtk DataFile Version 3.0\n')
f.write('Polydata Polygons\n')
f.write('ASCII\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS %d float\n'%(N**3))
np.savetxt(f, P, fmt = '%.3f')
f.write('POLYGONS %d %d\n'%(nF,5*nF))
np.savetxt(f, np.column_stack([4*np.ones(nF), F]), fmt='%d')

f.write('POINT_DATA %d\n'%(N**3))
f.write('VECTORS field float\n')
#f.write('LOOKUP_TABLE default\n')
np.savetxt(f, field, fmt = '%.3f')

f.close()
del f
#can generate streamlines
f = open('STRUCTURED_GRID.vtk','w')

f.write('# vtk DataFile Version 3.0\n')
f.write('Polydata Polygons\n')
f.write('ASCII\n')
f.write('DATASET STRUCTURED_GRID\n')
f.write('DIMENSIONS %d %d %d\n'%(N,N,N))
f.write('POINTS %d float\n'%(N**3))
np.savetxt(f, P, fmt = '%.3f')

f.write('POINT_DATA %d\n'%(N**3))
f.write('VECTORS field float\n')
np.savetxt(f, field, fmt = '%.3f')

f.close()
del f
########################################################################
########################################################################

