# cut_plane

import sys
import os
import numpy as np
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')
from scipy.integrate import odeint
from matplotlib.patches import Circle, PathPatch, Rectangle
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
import _CCMC as ccmc
import pos_sun as ps

#parameter_unit = 'nPa'
# Inputs for file and magnetic field model
year = 2003
day = 20
month = 11
hours = 7.
minutes = 0.
seconds = 0.

# Plot title
title = 'SCARR5 ' + str(year) + '-' + str(month) + '-' + str(day) + 'T07:00'
filename = conf["f_path"] + '3d__var_3_e' + str(year) + str(month) + str(day) + '-070000-000.out.cdf'

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

# run parameters
debug = False
Nb = 10
sign=-1  # changes sign of magnetic field used to trace the field lines

UT=hours*hr + minutes*minn + seconds*s
# Start point of main field line
MLON = 68.50*deg
MLAT = 50.00*deg


# open kameleon
kameleon = ccmc.Kameleon()
kameleon.open(filename)
print(filename, "Opened " + filename)
interpolator = kameleon.createNewInterpolator()

def ex_data(variable, x,y,z):
    # Get data from file, interpolate to point
    kameleon.loadVariable(variable)
    data = interpolator.interpolate(variable, x, y, z)
    return data

def data_in_U(variable,u,v,U1,U2):
    # Get the data in the U coordinates (defined by the cut plane vectors U1 and U2)
    x,y,z = u*U1+v*U2
    B=np.array([ex_data('bx', x,y,z), ex_data('by', x,y,z), ex_data('bz', x,y,z)])
    if variable == 'bu1':
        return np.dot(B,U1)
    if variable == 'bu2':
        return np.dot(B,U2)
    if variable == 'bu3':
        return np.dot(B,U3)
    else:
        return ex_data(variable, x,y,z)

def dXds(X, s):
    '''
    derivative function for field line ODE
    dx/ds = Bx(x,y,z)/Bm
    dy/ds = By(x,y,z)/Bm
    dz/ds = Bz(x,y,z)/Bm
    
    X = [x, y, z]
    B = [Bx, By, Bz]
    Bm = sqrt(Bx**2 + By**2 + Bz**2)
    s=arclength    
    '''
    B=np.array([ex_data('bx', X[0],X[1],X[2]), ex_data('by', X[0],X[1],X[2]), ex_data('bz', X[0],X[1],X[2])])
    Bm=np.sqrt(np.dot(B,B))
    if 1e-6<Bm<1e+6:
        return (sign/Bm)*B
    else:
        return [0., 0., 0.] #or nan

# main field line start point (in MAG)
phiMAG = MLON
thetaMAG = np.pi/2. - MLAT
R = 1.

u_st = R*np.sin(thetaMAG)*np.sin(phiMAG)
v_st = R*np.sin(thetaMAG)*np.cos(phiMAG)
w_st = R*np.cos(thetaMAG)
v = ps.MAGtoGSM([u_st, v_st, w_st], month, day, year, UT)
x_st = v[0]
y_st = v[1]
z_st = v[2]

if debug:
    print("---------")
    print(u_st)
    print(v_st)
    print(w_st)
    print("--------")
    print(x_st)
    print(y_st)
    print(z_st)

# Trace field lines
s_grid = np.linspace(0, 10., 100.)
solns = (np.nan)*np.empty((s_grid.size, 3, x_st.size))
X0 = [x_st, y_st, z_st] # Initial condition
soln = odeint(dXds, X0, s_grid)

# initialize vectors for defining field line cut plane
v1=(np.nan)*np.empty((3,))
v2=(np.nan)*np.empty((3,))
v3=(np.nan)*np.empty((3,))
U1=(np.nan)*np.empty((3,))
U2=(np.nan)*np.empty((3,))
U3=(np.nan)*np.empty((3,))
Mdipole=(np.nan)*np.empty((3,))

# define condition on the field line points
tr = np.logical_and(soln[:,0]**2+soln[:,1]**2+soln[:,2]**2 >=1., soln[:,0]**2+soln[:,1]**2+soln[:,2]**2 < 20.)
# create the arrays of the restricted field line componentwise
solx=soln[:,0]
solx=solx[tr]
soly=soln[:,1]
soly=soly[tr]
solz=soln[:,2]
solz=solz[tr]
# reasemble
sol=np.column_stack([solx,soly,solz])

# define vects for plain of main field line
v1 = sol[0,:]
v2 = sol[-1,:]
half = int(sol.shape[0]/2)
v3 = sol[half,:]
# define cut plane coordinates based on main field line 
# (U3 normal to the plane)
U2 = (v1-v2)/np.linalg.norm(v1-v2)
Mdipole = ps.MAGtoGSM([0.,0.,1.], month, day, year, UT)
U3 = np.cross(v3-v1, U2)/np.linalg.norm(np.cross(v3-v1, U2))
U1 = np.cross(U2, U3)
#return U1,U2,U3,Mdipole

#------------------------------
kameleon.close()
print("Closed " + filename)
#-------------------------------

f = open(conf["m_path"] + 'magnetosphere/data/cut_plane_info.txt','w')
f.write('%.7e %.7e %.7e\n'%(Mdipole[0], Mdipole[1], Mdipole[2]))
f.write('%.7e %.7e %.7e\n'%(U1[0], U1[1], U1[2]))
f.write('%.7e %.7e %.7e\n'%(U2[0], U2[1], U2[2]))
f.write('%.7e %.7e %.7e\n'%(U3[0], U3[1], U3[2]))

print (month,day,year,UT)
print 'Mdipole = ', Mdipole
print 'U1 = ', U1
print 'U2 = ', U2
print 'U3 = ', U3
