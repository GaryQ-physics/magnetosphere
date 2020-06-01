# sphere_vtk_write

import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

fname = '/home/gary/magnetosphere/vtk/rotated_sphere.vtk'
var='longitude'
Nt=300
Np=300

r2=(1/np.sqrt(2.))
Rot=np.array([[1, 0, 0], [0,r2, r2], [0,-r2, r2]])

print(Rot)

def func(x,y,z):
    return np.arctan2(y,x)

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

R=2.
theta = np.linspace(0.,np.pi,Nt)
phi = np.linspace(0.,2.*np.pi,Np)

#B2, B3, B1 =np.meshgrid(Y,Z,X)
#B1, B2, B3 =np.meshgrid(X,Y,Z)
#x=R*np.cos(phi)*np.sin(theta)
#y=R*np.sin(phi)*np.sin(theta)
#z=R*np.cos(theta)

B1, B2 = np.meshgrid(theta,phi)
B1= B1.flatten(order='C')
B2= B2.flatten(order='C')
B=np.column_stack((B1,B2))  # B[i,0]=B1[i], B[i,1]=B2[i]

x=R*np.cos(B[:,1])*np.sin(B[:,0])
y=R*np.sin(B[:,1])*np.sin(B[:,0])
z=R*np.cos(B[:,0])
XYZ=np.column_stack((x,y,z))


A=(np.nan)*np.empty((B1.size,))
xr=(np.nan)*np.empty((B1.size,))
yr=(np.nan)*np.empty((B1.size,))
zr=(np.nan)*np.empty((B1.size,))

for l in range(A.size):
    xr[l],yr[l],zr[l] = Rot.dot(XYZ[l,:])
    A[l]=func(XYZ[l,0], XYZ[l,1], XYZ[l,2])

XYZr=np.column_stack((xr,yr,zr))

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
f.close()
print("Wrote " + fname)
