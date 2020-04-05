# kameleon_field_paraview

#-----INSTRUCTIONS---------------------------
# run paraview
# File -> Open  kameleon_field_paraview.vtk
# click Apply
# under display(geometry representation), select points
# Filters -> Alphabetical -> Glyph
# click Apply
#--------------------------------------------

# Directory with kameleon subdirectory
#k_path = '/Users/robertweigel/'
k_path = '/home/gary/magnetosphere/'

# Directory of .cdf file
f_path = k_path + 'events/'

import sys
import numpy as np
# !!!path append needs to be ahead of import odeint for my computor!!!
sys.path.append(k_path + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(k_path + 'kameleon/lib/python2.7/site-packages/ccmc/')
from scipy.integrate import odeint
import _CCMC as ccmc

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

# Inputs for file and magnetic field model
year = 2003
day = 20
month = 11
hours = 7.
minutes = 0.
seconds = 0.

filename = f_path + '3d__var_3_e' + str(year) + str(month) + str(day) + '-070000-000.out.cdf'

# open kameleon
kameleon = ccmc.Kameleon()
kameleon.open(filename)
print(filename, "Opened " + filename)
interpolator = kameleon.createNewInterpolator()

def ex_data(variable, x,y,z):
    if np.sqrt(x**2+y**2+z**2)<1e-4: return 0.
    # Get data from file, interpolate to point
    kameleon.loadVariable(variable)
    data = interpolator.interpolate(variable, x, y, z)
    return data

N=10
fname = 'kameleon_field_paraview.vtk';

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
            f.write('%e %e %e\n'%( ex_data('bx',X[i],Y[j],Z[k]), ex_data('by',X[i],Y[j],Z[k]), ex_data('bz',X[i],Y[j],Z[k]) ))

f.close()

#------------------------------
kameleon.close()
print("Closed " + filename)
#-------------------------------
