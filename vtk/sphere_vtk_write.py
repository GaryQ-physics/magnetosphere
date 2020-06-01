# sphere_vtk_write

import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

fname = '/home/gary/magnetosphere/vtk/rotated_sphere.vtk'
var='longitude'
Nt = 100
Np = 100

r2=(1/np.sqrt(2.))
#Rot=np.array([[1, 0, 0], [0,r2, r2], [0,-r2, r2]])
Rot=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

print(Rot)

def func(x,y,z):
    return np.arctan2(y,x)

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

R = 1.
theta = np.linspace(0., np.pi, Nt)
phi = np.linspace(0., 2.*np.pi, Np)

B1, B2 = np.meshgrid(phi,theta)
B1 = B1.flatten(order='C')
B2 = B2.flatten(order='C')
B = np.column_stack((B1,B2))  # B[i,0]=B1[i], B[i,1]=B2[i]

normPhi = np.linspace(0., 1., Np)
normTheta = np.flipud(np.linspace(0., 1., Nt))
u, v = np.meshgrid(normPhi,normTheta)
u = u.flatten(order='C')
v = v.flatten(order='C')
UV = np.column_stack((u,v))

PI = np.pi*np.ones((B1.size, ))
x = R*np.cos(B1+PI)*np.sin(B2)
y = R*np.sin(B1+PI)*np.sin(B2)
z = R*np.cos(B2)
XYZ = np.column_stack((x, y, z))


A = (np.nan)*np.empty((B1.size, ))
xr = (np.nan)*np.empty((B1.size, ))
yr = (np.nan)*np.empty((B1.size, ))
zr = (np.nan)*np.empty((B1.size, ))

for l in range(A.size):
    xr[l],yr[l],zr[l] = Rot.dot(XYZ[l,:])
    A[l] = func(XYZ[l,0], XYZ[l,1], XYZ[l,2])

XYZr = np.column_stack((xr, yr, zr))

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Structured Grid for rotated sphere\n')
f.write('ASCII\n')
f.write('DATASET STRUCTURED_GRID\n')
f.write('DIMENSIONS ' + str(Nt) + ' ' + str(Np) + ' ' + str(1) + '\n' )
f.write('POINTS '+str(Nt*Np)+' float\n')
np.savetxt(f, XYZr)
f.write('\n')

f.write('POINT_DATA ' + str(Nt*Np) + '\n')
f.write('TEXTURE_COORDINATES TextureCoordinates 2 float\n')
np.savetxt(f, UV)

f.close()
print("Wrote " + fname)


