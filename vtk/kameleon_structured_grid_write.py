# kameleon_structured_grid_write

import sys
import os
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')
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

var = 'jy'
Ny = 3*30
Nz = 3*30
Nx_main = 3*27
Nx_tail = 30

filename = conf["f_path"] + '3d__var_3_e' + str(year) + str(month) + str(day) + '-070000-000.out.cdf'
fname = conf["m_path"] + 'magnetosphere/data/kameleon_structured_grid.vtk'

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

X = np.concatenate((np.linspace(-200,-20.05,Nx_tail), np.linspace(-20.,15.,Nx_main) ))
Nx = X.size
Y = np.linspace(-10.,10.,Ny)
Z = np.linspace(-10.,10.,Nz)
print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Structured Grid ' + var + '\n')
f.write('ASCII\n')
f.write('DATASET STRUCTURED_GRID\n')
f.write('DIMENSIONS ' + str(Nx) + ' ' + str(Ny) + ' ' + str(Nz) + '\n' )
f.write('POINTS '+str(Nx*Ny*Nz)+' float\n')
for k in range(Nz):
    for j in range(Ny):
        for i in range(Nx):
            f.write('%.2e %.2e %.2e\n'%(X[i], Y[j], Z[k]))
f.write('\n')
f.write('POINT_DATA ' + str(Nx*Ny*Nz) + '\n')
f.write('SCALARS point_scalars float 1\n')
f.write('LOOKUP_TABLE default\n')

for k in range(Nz):
    for j in range(Ny):
        for i in range(Nx):
            f.write('%.2e\n'%(ex_data(var, X[i], Y[j], Z[k])))
f.close()
print("Wrote " + fname)

#------------------------------
kameleon.close()
print("Closed " + filename)
#-------------------------------
