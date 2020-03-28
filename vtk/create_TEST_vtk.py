import numpy as np
#import os
#import sys
import matplotlib.pyplot as pl
#dname = '/home/wesley/Documents/CDS490_projects/cds-490/tsyganenko/'
#dname = os.pwd()

num_radii = 20

#Stepsize
step = .5

#Bounds of simulation
x_range = [-num_radii,num_radii,step]
y_range = [-num_radii,num_radii,step]
z_range = [-num_radii,num_radii,step]

#Name of output file
fname = 'dipole_VTK_test_' + str(num_radii) + '.vtk'

#Write command to file
print("Writing input_T01_01c.txt")
f = open('input_T01_01c.txt','w')
cmd = ('%.02f\n%.02f\n%.02f\n%.02f\n%.02f\n%.02f\n%.02f\n%.02f\n%.02f'%
(x_range[0],x_range[1],x_range[2],y_range[0],y_range[1],y_range[2],
 z_range[0],z_range[1],z_range[2]))
f.write(cmd)
f.close()

#Read command in, write Tsyganenko field data to file
#print("Writing T01_01c.txt")    
#com = './T01_01c < input_T01_01c.txt > T01_01c.txt'
#print("Executing " + com)
#os.system(com)

#Write input file
print("Writing input_dipole.txt")
f = open(r'input_dipole.txt','w')

f.write(cmd)
f.close()

#Read input file, write dipole magnetic field data to file.
#print("Writing dipole.txt")   
#com2 = './dipole < input_dipole.txt > dipole.txt'
#print("Executing " + com2)
#os.system(com2)

# Read both Tsyganenko and dipole data as variables.
#print("Reading T01_01c.txt")
#A = np.genfromtxt('T01_01c.txt',dtype=np.float)
#A[np.isnan(A)] = 0
#print("Reading dipole.txt")
#B = np.genfromtxt('dipole.txt',dtype=np.float)
#B[np.isnan(B)] = 0
A=np.zeros((81**3,4))
N = np.shape(A)[0]
num_divisions = np.ones(3)*np.shape((np.arange(-num_radii,num_radii+step,step)))
ARR=np.zeros((81**3,4))
B=np.zeros((81**3,4))
X=np.arange(x_range[0],x_range[1]+step,x_range[2])
Y=np.arange(y_range[0],y_range[1]+step,y_range[2])
Z=np.arange(z_range[0],z_range[1]+step,z_range[2])
l=0
R=0.
M = -31000.
for k in range(Z.size):
    for j in range(Y.size):
        for i in range(X.size):
            R=np.sqrt(X[i]**2+Y[j]**2+Z[k]**2)
            if(R>0.01):
                B[l,0]=M*(3*X[i]**2 - R**2)/R**5
                B[l,1]=3*M*Y[j]*X[i]/R**5
                B[l,2]=3*M*Z[k]*X[i]/R**5
                B[l,3]=R
            else:
                B[l,0]=0.
                B[l,1]=0.
                B[l,2]=0.
                B[l,3]=0.
            l=l+1

ARR[0,0]=0.





#Fake density values with r-dependency only.
#density = 5e11/B[:,3]
#density[B[:,3] < 1.0] = 5e11
density=5e11*np.ones(B.shape[0],)

#Neutral sheet definition taken from "Regions and Boundaries in the SSC Software:
#Algorithms"
sheet = np.zeros((N*2,3))
psi = np.pi/3.
Rh=8.
d=4.
G=10.
Ly=10.
#Calculate neutral sheet.
index=0
#for x in np.arange(-num_radii,num_radii+1,x_range[2]):
#    for y in np.arange(-num_radii,num_radii+1,y_range[2]):
#        z1 = (0.5*np.tan(psi) * np.sqrt((x-Rh*np.cos(psi))**2 + (d*np.cos(psi))**2)
#            - np.sqrt((x+Rh*np.cos(psi))**2 + (d*np.cos(psi))**2)
#            - G*np.sin(psi) * y**4/(y**4+Ly**4))
#        z2 = -z1
#        print str([x,y,z1])
#        sheet[index,:] = [x,0,z1]
#        index += 1
#        sheet[index,:] = [x,0,z2]
#        index += 1
for x in np.arange(-num_radii,num_radii+1,x_range[2]):
    z1 = 0.5*np.tan(psi) * (np.sqrt((x-Rh*np.cos(psi))**2 + (d*np.cos(psi))**2) 
        - np.sqrt((x+Rh*np.cos(psi))**2 + (d*np.cos(psi))**2))
    z2 = -z1
    sheet[index,:] = [x,0,z1]
    index += 1
    sheet[index,:] = [x,0,z2]
    index += 1

print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('T01_01 and Dipole vectors\n')
f.write('ASCII\n')
f.write('\n')
f.write('DATASET STRUCTURED_POINTS\n')
f.write('DIMENSIONS %i %i %i\n'%(num_divisions[0],num_divisions[1],num_divisions[2]))
f.write('ORIGIN %f %f %f\n'%(-num_radii,-num_radii,-num_radii))
f.write('SPACING %f %f %f\n'%(x_range[2],y_range[2],z_range[2]))

f.write('\n')
f.write('POINT_DATA %i\n'%N)
f.write('VECTORS T01_01 float\n')

for i in range(N):
    f.write('%e %e %e\n'%(A[i,0],A[i,1],A[i,2]))
f.write('\n')
f.write('VECTORS DIPOLE float\n')
for i in range(N):
    f.write('%e %e %e\n'%(B[i,0],B[i,1],B[i,2]))
f.write('\n')
f.write('SCALARS mock_ionosphere float\n')
f.write('LOOKUP_TABLE default\n')
    
for i in range(N):
    f.write('%e\n'%(density[i]))
f.close()

#Produce neutral sheet
point_num = sheet.shape[0]

f = open('neutralsheet.vtk','w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Neutral sheet\n')
f.write('ASCII\n')
f.write('\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS %i float\n'%point_num)
for i in range(point_num):
    f.write('%0.02f %0.02f %0.02f\n'%(sheet[i,0],sheet[i,1],sheet[i,2]))
f.close()
pl.plot(sheet[:,0],sheet[:,2],'r.')
