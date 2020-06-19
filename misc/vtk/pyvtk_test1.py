"""
pip install pyvtk
"""

# pyvtk_test1

#import sys
#sys.path = ['..']+sys.path
import numpy as np
import pyvtk
import time

npscal = True
npgrid = False
# fastest appears to be npscal = True , npgrid = False
# although not much difference


a = 60
b = 80
c = 120

# make grid (2 formats)------
X = np.array([i for i in range(a)])
Y = np.array([j for j in range(b)])
Z = np.array([k for k in range(c)])
B2, B3, B1 = np.meshgrid(Y, Z, X)
#B1, B2, B3 = np.meshgrid(X, Y, Z) # seems more natural but doesnt work with vtk structured_grid format
B1 = B1.flatten(order='C')
B2 = B2.flatten(order='C')
B3 = B3.flatten(order='C')

p_np = np.column_stack((B1, B2, B3))
p_tup = [(i,j,k) for k in range(c) for j in range(b) for i in range(a)]
#------------------


# def scalars by function and then make np array ------------------------
def f(x,y,z):
    return x*y*z

fvals = np.nan*np.empty((p_np.shape[0],))
for k in range(p_np.shape[0]):
    fvals[k]=f(p_np[k,0], p_np[k,1], p_np[k,2])
#----------------------

to = time.time()

if npscal:
    pointdata = pyvtk.PointData(\
    pyvtk.Scalars(fvals,
            name='sample_scalars',
            lookup_table='my_table'),
            )
    if npgrid:
        vtk = pyvtk.VtkData(pyvtk.StructuredGrid([a,b,c],p_np), pointdata)
    else:
        vtk = pyvtk.VtkData(pyvtk.StructuredGrid([a,b,c],p_tup), pointdata)
    
else:
    if npgrid:
        vtk = pyvtk.VtkData(pyvtk.StructuredGrid([a,b,c],p_np))
    else:
        vtk = pyvtk.VtkData(pyvtk.StructuredGrid([a,b,c],p_tup))
    vtk.point_data.append(vtk.structure.Scalars(f,'x*y*z'))

vtk.tofile('pyvtk_test1','ascii')

tf = time.time()

print('--------------')
print('time = ', tf-to)
print('npscal = ', npscal)
print('npgrid = ', npgrid)
print('--------------')
