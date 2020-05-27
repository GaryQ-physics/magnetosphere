# cut_plane_python_plotting

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config_paths import config
conf = config()

import matplotlib.pyplot as plt
# Following needed for projection='3d', but PyFlakes will warn that not used
from scipy.integrate import odeint

sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')
import _CCMC as ccmc
import pos_sun as ps

from cut_plane import ex_data, dXds

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

# run parameters
debug = False
n=50 # number of pts on cutplane grid 
m=50
# parameter to plot
parameter =  'p'

# Start point of main field line
mlon = 68.50
mlat = 50.00
r = 1.01

# Inputs for file and magnetic field model
year = 2003
day = 20
month = 11
hours = 7
minutes = 0
seconds = 0


#def plot(Event, parameter, xlim=xlim, ylim=ylim, nx=n, ny=m, png=True):

def data_in_U(kam, interp, variable, u, v, U1, U2, U3):
    # Get the data in the U coordinates (defined by the cut plane vectors U1 and U2)
    x,y,z = u*U1+v*U2
    B = np.array([ex_data(kam,interp,'bx', x,y,z), 
                  ex_data(kam,interp,'by', x,y,z), 
                  ex_data(kam,interp,'bz', x,y,z)])
    if variable == 'bu1':
        return np.dot(B,U1)
    if variable == 'bu2':
        return np.dot(B,U2)
    if variable == 'bu3':
        return np.dot(B,U3)
    else:
        return ex_data(kam,interp, variable, x,y,z)


time = [year, month, day, hours, minutes, seconds]

# Plot title
title = 'SCARR5 ' + str(year) + '-' + str(month) + '-' + str(day) + 'T07:00'
title = title + "\n" + "[mlat,mlon]=[{0:.1f}, {1:.1f}]".format(mlat, mlon)

filename = conf["run_path"] + '3d__var_3_e' \
            + '%04d%02d%02d-%02d%02d%02d-000' % tuple(time)

filename_in = filename + '.out.cdf'
filename_out = filename + '.png'

X0 = ps.MAGtoGSM([r, mlat, mlon], time, 'sph', 'car')

# open kameleon
kameleon = ccmc.Kameleon()
print("Opening " + filename_in)
kameleon.open(filename_in)
print("Opened " + filename_in)
interpolator = kameleon.createNewInterpolator()

parameter_unit = kameleon.getVisUnit(parameter)

# Trace field line
s_grid = np.linspace(0., 10., 100)
soln = odeint(dXds, X0, s_grid, args=(kameleon, interpolator))
if debug:
    print(X0)
    print(np.dot(X0,X0))
    print(soln)

# initialize vectors for defining field line cut plane
v1=(np.nan)*np.empty((3,))
v2=(np.nan)*np.empty((3,))
v3=(np.nan)*np.empty((3,))
U1=(np.nan)*np.empty((3,))
U2=(np.nan)*np.empty((3,))
U3=(np.nan)*np.empty((3,))
Mdipole=(np.nan)*np.empty((3,))

# define restriction condition on the field line points
tr = np.logical_and(soln[:,0]**2 + soln[:,1]**2 + soln[:,2]**2 >= 1.,
                    soln[:,0]**2 + soln[:,1]**2 + soln[:,2]**2 < 20.)

# restricted field lines
sol = soln[tr,:]

# define vects for plane of main field line
v1 = sol[0,:]  # First point on field line
v2 = sol[-1,:] # Last point on field line
half = int(sol.shape[0]/2) 
v3 = sol[half,:] # Approximate mid-point on field line

# define cut plane coordinates based on main field line 
# (U3 is normal to the plane)
U2 = (v1-v2)/np.linalg.norm(v1-v2)
Mdipole = ps.MAGtoGSM([0., 0., 1.], time, 'car', 'car')
U3 = np.cross(v3-v1, U2)/np.linalg.norm(np.cross(v3-v1, U2))
U1 = np.cross(U2, U3)

x_1d = np.linspace(0, 4, n)
y_1d = np.linspace(-3, 3, m)
X, Y = np.meshgrid(x_1d, y_1d) # grid of points on the cutplane
Z = np.zeros((n, m))
for i in range(n):
    for j in range(m):
        # grid of the corresponding values of variable. To be color plotted
        Z[i,j] = data_in_U(kameleon, interpolator, parameter, X[i,j], Y[i,j], U1, U2, U3)

kameleon.close()
print("Closed " + filename_in + "\n")
print('Close plot to exit')

# Plotting
plt.clf()
fig = plt.figure(dpi=200)
ax2 = fig.add_subplot(1,2,2)

# Plot cut plane data
ax2.set(title=title)
ax2.set(xlabel="Tailward distance [$R_E$]")
ax2.set(ylabel="Northward distance [$R_E$]")
ax2.axis('square')
pcm = ax2.pcolormesh(X, Y, Z)

# Reason for choice of fraction and pad:
# https://stackoverflow.com/a/39948312/1491619
cb = fig.colorbar(pcm, ax=ax2, fraction=0.04, pad=0.08)
cb.ax.set_title(parameter + ' [' + parameter_unit + ']', ha='left')
plt.xlim(0, 4)
plt.ylim(-3, 3)

# Add field line to 2D plot        
solCut=np.zeros((sol.shape[0],2))
for k in range(sol.shape[0]):
    solCut[k,0]=np.dot(sol[k,:],U1)
    solCut[k,1]=np.dot(sol[k,:],U2)
ax2.plot(solCut[:,0],solCut[:,1], 'red', lw=1)


plt.show()
