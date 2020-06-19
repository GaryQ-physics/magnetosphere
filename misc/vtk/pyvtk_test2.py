# pyvtk_test2

#import sys
#sys.path = ['..']+sys.path
import numpy as np
import pyvtk
import time

npgrid = False


a = 6
b = 8
c = 12

def f(x,y,z):
    return x*y*z

to = time.time()

# -- ----------
if npgrid:
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
else:
    p_tup = [(i,j,k) for k in range(c) for j in range(b) for i in range(a)]
    fvals = [f(*p_tup[l]) for l in range(len(p_tup))] # https://stackoverflow.com/questions/32896651/pass-multiple-arguments-in-form-of-tuple
#------------------


pointdata = pyvtk.PointData(\
pyvtk.Scalars(fvals,
        name='sample_scalars',
        lookup_table='my_table'),
        )
if npgrid:
    vtk = pyvtk.VtkData(pyvtk.StructuredGrid([a,b,c],p_np), pointdata)
else:
    vtk = pyvtk.VtkData(pyvtk.StructuredGrid([a,b,c],p_tup), pointdata)



vtk.tofile('pyvtk_test2','ascii')
#vtk.tofile('/home/gary/magnetosphere/test_data/pyvtk_sg_test2','binary')


tf = time.time()

print('--------------')
print('num points is ', a*b*c)
print('time = ', tf-to)
print('npgrid = ', npgrid)
print('--------------')
