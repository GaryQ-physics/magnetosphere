# pyvista only works with python3
import numpy as np
import vtk
import pyvista as pv

N=10

xax = yax = zax = 1.*np.arange(N)
#x, y, z = np.meshgrid(ax,ax,ax)
x, y = np.meshgrid(xax,yax)
x = x.flatten()
y = y.flatten()
assert(np.all(x == np.kron(np.ones(N),xax)))
assert(np.all(y == np.kron(yax,np.ones(N))))

#xy = np.column_stack([x,y])

x == np.kron(np.ones(N),xax)
y == np.kron(yax,np.ones(N))

z = np.kron(zax, np.ones(N**2))

x = np.kron(np.ones(N), x)
y = np.kron(np.ones(N), y)

P = np.column_stack([x,y,z])

'''
E_temp = [(0,1),(1,2),
     (3,4),(4,5),
     (6,7),(7,8),
     (0,3),(1,4),(2,5),
     (3,6),(4,7),(5,8)]
E = []
for i in range(N-1):
    for j in range(N):
        E.append((ind[i,j], ind[i+1,j]))
for i in range(N):
    for j in range(N-1):
        E.append((ind[i,j], ind[i,j+1]))
assert(set(E)==set(E_temp))
'''

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

f.write('POINT_DATA %d\n'%(N**3))
f.write('VECTORS field float\n')
np.savetxt(f, field, fmt = '%.3f')

f.close()
del f

#pyvista outbuts "# vtk DataFile Version 5.1"  instead of 3.0 , and it has extra sections 
cells = np.column_stack([8*np.ones(nV, dtype=int), V]).flatten()
offset = (8+1)*np.arange(nV)
print(nV)
print(cells)
print(offset)
celltypes = np.empty(nV, dtype=np.uint8)
celltypes[:] = vtk.VTK_HEXAHEDRON

print(offset.dtype)
print(cells.dtype)
print(celltypes.dtype)
print(P.dtype)

grid = pv.UnstructuredGrid(offset, cells, celltypes, P)
grid.save('UNSTRUCTURED_GRID-volumes-pyvista.vtk', binary=False)


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
f.write('CELLS %d %d\n'%(nE,4*nE))
np.savetxt(f,np.column_stack([3*np.ones(nE), E]), fmt='%d')
f.write('CELL_TYPES %d\n'%(nE))
np.savetxt(f, 3*np.ones(nE), fmt = '%d')

f.write('POINT_DATA %d\n'%(N**3))
f.write('VECTORS field float\n')
np.savetxt(f, field, fmt = '%.3f')

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

f.write('POINT_DATA %d\n'%(N**3))
f.write('VECTORS field float\n')
np.savetxt(f, field, fmt = '%.3f')

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

