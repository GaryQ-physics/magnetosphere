# coordinate_test

# Directory with kameleon subdirectory
#k_path = '/Users/robertweigel/'
k_path = '/home/gary/magnetosphere/'

# Directory of .cdf file
f_path = './'

import sys
import numpy as np

# !!!path append needs to be ahead of import odeint for my computor!!!
sys.path.append(k_path + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(k_path + 'kameleon/lib/python2.7/site-packages/ccmc/')
sys.path.append(k_path + 'events/')
# For OS-X
import matplotlib as mpl
mpl.use('TkAgg')

import matplotlib.pyplot as plt
# Needed for projection='3d', but Spyder
# will warn that not used
from mpl_toolkits import mplot3d
from scipy.integrate import odeint
from matplotlib.patches import Circle, PathPatch, Rectangle
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
import time

import _CCMC as ccmc
import pos_sun as ps

# Inputs for file and magnetic field model
year = 2003
day = 20
month = 11
hours = 7.
minutes = 0.
seconds = 0.

# Plot title
title = 'SCARR5 ' + str(year) + '-' + str(month) + '-' + str(day) + 'T07:00'
filename = f_path + '3d__var_3_e' + str(year) + str(month) + str(day) + '-070000-000.out.cdf'

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

# run parameters
debug = False
usePatch = True
sign=-1  # changes sign of magnetic field used to trace the field lines
l=4
UT=hours*hr + minutes*minn + seconds*s
#epochTime = long(1069311600) # from https://www.epochconverter.com/



# open kameleon
kameleon = ccmc.Kameleon()
kameleon.open(filename)
print(filename, "Opened " + filename)
interpolator = kameleon.createNewInterpolator()
coordinate_interpolator = kameleon.createCoordinateInterpolator() # no arguments assumes native
print 'epoch time:', coordinate_interpolator.getEphemTime(), 'seconds'
coordinate_interpolator.setEphemTime(1069311600)
print 'epoch time:', coordinate_interpolator.getEphemTime(), 'seconds'
print coordinate_interpolator.get_model_coords()

def ex_data(variable, x,y,z):
    # Get data from file, interpolate to point
    kameleon.loadVariable(variable)
    data = interpolator.interpolate(variable, x, y, z)
    return data

R=3.
u = R*np.random.rand(l,)
v = R*np.random.rand(l,)
w = R*np.random.rand(l,)
x = np.nan*np.empty(l,)
y = np.nan*np.empty(l,)
z = np.nan*np.empty(l,)
print 'u=', u
print 'v=', v
print 'w=', w

print 'x=', x
print 'y=', y
print 'z=', z

print coordinate_interpolator.get_model_coords()
coordinate_interpolator.setPreferredCoordinates('MAG')
print coordinate_interpolator.get_model_coords()
for i in range(l):
    print i
    print 'l 1'
    MAG_coords = ccmc.Position()
    print 'l 2'
    MAG_coords.c0 = u[i]
    print 'l 3'
    MAG_coords.c1 = v[i]
    print 'l 4'
    MAG_coords.c2 = w[i]
    print 'l 5'
    GSM_coords_ieModel = ccmc.Position()
    print 'l 6'
    coordinate_interpolator.convertCoordinates(MAG_coords, GSM_coords_ieModel)
    print 'l 7'
    x[i] = GSM_coords_ieModel.c0
    y[i] = GSM_coords_ieModel.c1
    z[i] = GSM_coords_ieModel.c2


print 'x=', x
print 'y=', y
print 'z=', z

print'mark1'

#coordinate_interpolator.setPreferredCoordinates("GSM")
x=np.nan*x
y=np.nan*y
z=np.nan*z
print'mark2'
print 'u=', u
print 'v=', v
print 'w=', w

print 'x=', x
print 'y=', y
print 'z=', z

for i in range(l):
    out = ps.MAGtoGSM([u[i], v[i], w[i]], month, day, year, UT)
    x[i] = out[0]
    y[i] = out[1]
    z[i] = out[2]
print 'mark3'
print 'x=', x
print 'y=', y
print 'z=', z
#print MAG_coords.c0, MAG_coords.c1, MAG_coords.c2
#print GSM_coords_ieModel.c0, GSM_coords_ieModel.c1, GSM_coords_ieModel.c2

#------------------------------
#kp.k_close()
kameleon.close()
print("Closed " + filename)
#-------------------------------
