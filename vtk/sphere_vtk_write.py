# sphere_vtk_write

import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

fname = conf["m_path"] + 'magnetosphere/data/sphere_vtk.vtk'
var='longitude'
Nt=300
Np=300

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
for l in range(A.size):
    A[l]=func(XYZ[l,0], XYZ[l,1], XYZ[l,2])

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Structured Grid for sphere with ' + var + '\n')
f.write('ASCII\n')
f.write('DATASET STRUCTURED_GRID\n')
f.write('DIMENSIONS ' + str(Nt) + ' ' + str(Np) + ' ' + str(1) + '\n' )
f.write('POINTS '+str(Nt*Np)+' float\n')
np.savetxt(f, XYZ)
f.write('\n')
f.write('POINT_DATA ' + str(Nt*Np) + '\n')
f.write('SCALARS point_scalars float 1\n')
f.write('LOOKUP_TABLE default\n')
np.savetxt(f, A)
f.close()
print("Wrote " + fname)
