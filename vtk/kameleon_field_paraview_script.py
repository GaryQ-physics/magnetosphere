# Viewing output file in Parview
# Open Paraview
# Select File -> Open -> kameleon_field_paraview.vtk
# Click Apply
# Under display(geometry representation), select points
# Filters -> Alphabetical -> Glyph
# Click Apply

import os
import sys
import numpy as np
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
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
N = 10

filename = conf['f_path'] + '3d__var_3_e' + str(year) + str(month) + str(day) + '-070000-000.out.cdf'
fname = '../data/kameleon_field_paraview.vtk';

if not os.path.exists(filename):
    print('Did not find ', filename)
    sys.exit(1)

# open kameleon
kameleon = ccmc.Kameleon()
try:
    kameleon.open(filename)
    print("Opened " + filename)
except:
    print('kameleon.open(' + filename + ') failed.')
    sys.exit(1)

interpolator = kameleon.createNewInterpolator()

def ex_data(variable, x,y,z):
    if np.sqrt(x**2+y**2+z**2)<1e-4: return 0.
    # Get data from file, interpolate to point
    kameleon.loadVariable(variable)
    data = interpolator.interpolate(variable, x, y, z)
    return data

X = np.linspace(-5.,5.,N)
Y = np.linspace(-5.,5.,N)
Z = np.linspace(-5.,5.,N)
print('length = ' + str(X.size))
print("Writing " + fname)
f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Data at unstructured points\n')
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
            f.write('%e %e %e\n' % 
                    (ex_data('bx',X[i],Y[j],Z[k]),
                     ex_data('by',X[i],Y[j],Z[k]),
                     ex_data('bz',X[i],Y[j],Z[k])))

f.close()
print("Wrote " + fname)
#------------------------------
kameleon.close()
print("Closed " + filename)
#-------------------------------
