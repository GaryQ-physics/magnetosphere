# write_field_test 

#-----INSTRUCTIONS---------------------------
# run paraview
# File -> Open  write_field_test.vtk
# click Apply
# under display(geometry representation), select points
# Filters -> Alphabetical -> Glyph
# click Apply
#--------------------------------------------

import numpy as np

def ex_data(variable, x,y,z):
    M=np.array([0.5,-0.5,1])  #dipole moment
    R=np.array([x,y,z])
    r=np.sqrt(np.dot(R,R))
    if r<1e-4:
        return 0.
    B=(1./(4*np.pi))*( (3./r**5)*np.dot(M,R)*R - (1./r**3)*M )
    data = B[variable]
    return data

N=10
fname = 'write_field_test.vtk';

X=np.linspace(-5.,5.,N)
Y=np.linspace(-5.,5.,N)
Z=np.linspace(-5.,5.,N)
print('length=',X.size)
print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Wind velocity at unstructured points\n')
f.write('ASCII\n')
f.write('\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS '+str(N**3)+' float\n')
for i in range(N):
    for j in range(N):
        for k in range(N):
            f.write('%e %e %e\n'%(X[i],Y[j],Z[k]))


f.write('\n')
f.write('POINT_DATA '+str(N**3)+'\n')
f.write('VECTORS point_vectors float\n')
for i in range(N):
    for j in range(N):
        for k in range(N):
            f.write('%e %e %e\n'%( ex_data(0,X[i],Y[j],Z[k]), ex_data(1,X[i],Y[j],Z[k]), ex_data(2,X[i],Y[j],Z[k]) ))


f.close()
