# cut_plane_python_plotting

import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config_paths import config
conf = config()
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')

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
import _CCMC as ccmc
import pos_sun as ps

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

# run parameters
debug = False
usePatch = False
plot3d = False
Nb = 10
sign=-1  # changes sign of magnetic field used to trace the field lines
n=50 # number of pts on cutplane grid 
m=50
# parameter to plot
parameter =  'p'
parameter_unit = 'nPa'

def ex_data(kam, interp, variable, x,y,z):
    # Get data from file, interpolate to point
    kam.loadVariable(variable)
    data = interp.interpolate(variable, x, y, z)
    if( x**2 + y**2 + z**2 >=1.):
        return data
    else:
        return 0.

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

def dXds(X, s, kam, interp):
    """
    derivative function for field line ODE
    dx/ds = Bx(x,y,z)/Bm
    dy/ds = By(x,y,z)/Bm
    dz/ds = Bz(x,y,z)/Bm
    
    X = [x, y, z]
    B = [Bx, By, Bz]
    Bm = sqrt(Bx**2 + By**2 + Bz**2)
    s=arclength    
    """
    
    B = np.array([ex_data(kam,interp, 'bx', X[0],X[1],X[2]), 
                  ex_data(kam,interp, 'by', X[0],X[1],X[2]), 
                  ex_data(kam,interp, 'bz', X[0],X[1],X[2])])
    Bm = np.sqrt(np.dot(B,B))
    if 1e-9 < Bm <1e+7:
        return (sign/Bm)*B
    else:
        if debug:
            if(Bm >= 1e+7): print('FIELD TOO HIGH')
            if(Bm <= 1e-7): print('FIELD TOO LOW')
        return [0., 0., 0.] # TODO: Return np.nan?

# functions from https://stackoverflow.com/questions/18228966/how-can-matplotlib-2d-patches-be-transformed-to-3d-with-arbitrary-normals
def rotation_matrix(d):
    """
    Calculates a rotation matrix given a vector d. The direction of d
    corresponds to the rotation axis. The length of d corresponds to 
    the sin of the angle of rotation.

    Variant of: http://mail.scipy.org/pipermail/numpy-discussion/2009-March/040806.html
    """
    sin_angle = np.linalg.norm(d)

    if sin_angle == 0:
        return np.identity(3)

    d /= sin_angle

    eye = np.eye(3)
    ddt = np.outer(d, d)
    skew = np.array([[    0,  d[2],  -d[1]],
                  [-d[2],     0,  d[0]],
                  [d[1], -d[0],    0]], dtype=np.float64)

    M = ddt + np.sqrt(1 - sin_angle**2) * (eye - ddt) + sin_angle * skew
    return M

def pathpatch_2d_to_3d(pathpatch, z = 0, normal = 'z'):
    """
    Transforms a 2D Patch to a 3D patch using the given normal vector.

    The patch is projected into they XY plane, rotated about the origin
    and finally translated by z.
    """
    if type(normal) is str: #Translate strings to normal vectors
        index = "xyz".index(normal)
        normal = np.roll((1.0,0,0), index)

    normal = (1./np.linalg.norm(normal))*normal #Make sure the vector is normalised

    path = pathpatch.get_path() #Get the path and the associated transform
    trans = pathpatch.get_patch_transform()

    path = trans.transform_path(path) #Apply the transform

    pathpatch.__class__ = art3d.PathPatch3D #Change the class
    pathpatch._code3d = path.codes #Copy the codes
    pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color    

    verts = path.vertices #Get the vertices in 2D

    d = np.cross(normal, (0, 0, 1)) #Obtain the rotation vector   
    M = rotation_matrix(d) #Get the rotation matrix

    pathpatch._segment3d = np.array([np.dot(M, (x, y, 0)) + (0, 0, z) for x, y in verts])

def pathpatch_translate(pathpatch, delta):
    """
    Translates the 3D pathpatch by the amount delta.
    """
    pathpatch._segment3d += delta

# Inputs for file and magnetic field model
year = 2003
day = 20
month = 11
hours = 7
minutes = 0
seconds = 0

time = [year,month,day,hours,minutes,seconds]
T=tuple(time)

# Plot title
title = 'SCARR5 ' + str(year) + '-' + str(month) + '-' + str(day) + 'T07:00'
filename = conf["run_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % T + '.out.cdf'

UT=hours*hr + minutes*minn + seconds*s

# Start point of main field line
MLON = 68.50
MLAT = 50.00

R = 1.01
X0 = ps.MAGtoGSM([R,MLAT,MLON],time,'sph','car')

# open kameleon
kameleon = ccmc.Kameleon()
kameleon.open(filename)
print(filename, "Opened " + filename)
interpolator = kameleon.createNewInterpolator()


# Trace field line
s_grid = np.linspace(0., 10., 100.)
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

# define condition on the field line points
tr = np.logical_and(soln[:,0]**2+soln[:,1]**2+soln[:,2]**2 >=1., soln[:,0]**2+soln[:,1]**2+soln[:,2]**2 < 20.)

# create the arrays of the restricted field line componentwise
solx=soln[:,0]
solx=solx[tr]
soly=soln[:,1]
soly=soly[tr]
solz=soln[:,2]
solz=solz[tr]

# reassemble
sol=np.column_stack([solx,soly,solz])
#print(sol)

# define vects for plane of main field line
v1 = sol[0,:]
v2 = sol[-1,:]
half = int(sol.shape[0]/2)
v3 = sol[half,:]

# define cut plane coordinates based on main field line 
# (U3 is normal to the plane)
U2 = (v1-v2)/np.linalg.norm(v1-v2)
Mdipole = ps.MAGtoGSM([0.,0.,1.],time,'car','car')
U3 = np.cross(v3-v1, U2)/np.linalg.norm(np.cross(v3-v1, U2))
U1 = np.cross(U2, U3)


x_1d = np.linspace(0, 4, n)
y_1d = np.linspace(-3, 3, m)
X, Y = np.meshgrid(x_1d, y_1d) # grid of points on the cutplane
Z = np.zeros((n, m))
for i in range(n):
    for j in range(m):
        # grid of the corresponding values of variable. To be color plotted
        Z[i,j]=data_in_U(kameleon, interpolator, parameter, X[i,j], Y[i,j], U1, U2, U3)

#------------------------------
#kp.k_close()
kameleon.close()
print("Closed " + filename)
#-------------------------------


print('Close plot to exit')

# Plotting
plt.clf()
fig = plt.figure(figsize=plt.figaspect(0.5))
if plot3d: ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax2 = fig.add_subplot(1,2,2)
#ax = fig.add_subplot(111, projection='3d')

# Plot cut plane data in right panel
ax2.set(title=title)
ax2.set(xlabel="Tailward distance [$R_E$]")
ax2.set(ylabel="Northward distance [$R_E$]")
ax2.axis('square')
pcm = ax2.pcolormesh(X, Y, Z)

# Reason for choice of fraction and pad:
# https://stackoverflow.com/a/39948312/1491619
cb = fig.colorbar(pcm, ax=ax2, fraction=0.046, pad=0.04)
cb.ax.set_title(parameter + ' [' + parameter_unit + ']')
#cb.ax.set_ylabel('nPa')
plt.xlim(0, 4)
plt.ylim(-3, 3)

# Plot field lines
if plot3d:
    for i in range(Nb+1): 
        from_list=solns_restr[i]
        sol=np.array(from_list)
        if (i == 0):
            # Event field line
            solCut=np.zeros((sol.shape[0],2))
            for k in range(sol.shape[0]):
                solCut[k,0]=np.dot(sol[k,:],U1)
                solCut[k,1]=np.dot(sol[k,:],U2)
            # Add field line to 3D plot
            ax1.plot3D(sol[:,0], sol[:,1], sol[:,2], 'red', lw=4)
            # Add field line to 2D plot        
            ax2.plot(solCut[:,0],solCut[:,1], 'red', lw=1)
        else:
            ax1.plot3D(sol[:,0], sol[:,1], sol[:,2], 'gray')

    ax1.view_init(elev=-120., azim=-16)
    ax1.set(xlabel="$X/R_E$ (GSM)")
    ax1.set(ylabel="$Y/R_E$ (GSM)")
    ax1.set(zlabel="$Z/R_E$ (GSM)")
    ax1.axis('square')

    ax1.set_xlim(-4, 4)
    ax1.set_ylim(-4, 4)
    ax1.set_zlim(-4, 4)

if plot3d:
    # Plot plane of field line and x-z plane asociated with y location of event.
    # define the cut plane (as span U1, U2) and the plane parallel to x-z plane
    para1, para2 = np.meshgrid(np.linspace(-1., 3., n), np.linspace(-2., 2., m))
    X_slice = U1[0]*para1+U2[0]*para2+v1[0]*np.ones((n,m))
    Y_slice = U1[1]*para1+U2[1]*para2+v1[1]*np.ones((n,m))
    Z_slice = U1[2]*para1+U2[2]*para2+v1[2]*np.ones((n,m))

    X_o = -1*para1 + 0*para2 + v1[0]*np.ones((n, m))
    Y_o =  0*para1 + 0*para2 + v1[1]*np.ones((n, m))
    Z_o =  0*para1 + 1*para2 + v1[2]*np.ones((n, m))

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi/2., 100)
    z_day = np.outer(np.cos(u), np.sin(v))
    y_day = np.outer(np.sin(u), np.sin(v))
    x_day = np.outer(np.ones(np.size(u)), np.cos(v))
    v = np.linspace(np.pi/2., np.pi, 100)
    z_night = np.outer(np.cos(u), np.sin(v))
    y_night = np.outer(np.sin(u), np.sin(v))
    x_night = np.outer(np.ones(np.size(u)), np.cos(v))

    #plot the planes in the 3D plot
    from matplotlib import cm as cmap
    if usePatch == False:
        ax1.plot_surface(X_slice, Y_slice, Z_slice, color='g', rstride=1, cstride=1,
                                linewidth=0, antialiased=False)
        ax1.plot_surface(X_o, Y_o, Z_o, color='b', rstride=1, cstride=1,
                                linewidth=0, antialiased=False)
    ax1.plot_surface(x_day, y_day, z_day, color='w', rstride=1, cstride=1,
                           linewidth=0, antialiased=False)
    ax1.plot_surface(x_night, y_night, z_night, color='black', rstride=1, cstride=1,
                           linewidth=0, antialiased=False)
    #ax.set(xlabel='x', ylabel='y', zlabel='z')
    #p = Circle((1, 1), 1)



    if usePatch:
        p = Rectangle((-2, -4), 4, 8, color='b')
        q = Rectangle((-2, -4), 4, 8, color='g')
        ax1.add_patch(p)
        ax1.add_patch(q)
        n=np.array([0,1,0])
        m=U3
        pathpatch_2d_to_3d(p, 0., n)
        pathpatch_2d_to_3d(q, 0., m)
        pathpatch_translate(p, v1)
        pathpatch_translate(q, v1)

plt.show()
