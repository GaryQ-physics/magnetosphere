# field_line_vtk_script



# Directory with kameleon subdirectory
#k_path = '/Users/robertweigel/'
k_path = '/home/gary/magnetosphere/'

# Directory of .cdf file
f_path = k_path + 'events/'

import sys
import numpy as np
sys.path.append(k_path + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(k_path + 'kameleon/lib/python2.7/site-packages/ccmc/')
sys.path.append(k_path + 'events/')
from scipy.integrate import odeint
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
Nb = 10
sign=-1  # changes sign of magnetic field used to trace the field lines
n=50 # number of pts on cutplane grid 
m=50

UT=hours*hr + minutes*minn + seconds*s

parameter='p'
# Start point of main field line
MLON = 68.50*deg
MLAT = 50.00*deg
fname = 'field_line.vtk'


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

# Background field lines start points (in MAG)
phi_st = np.linspace(0, 2*np.pi, Nb)
theta = np.pi/2.
R = 3.
u_st = R*np.sin(theta)*np.sin(phi_st)
v_st = R*np.sin(theta)*np.cos(phi_st)
w_st = R*np.cos(theta)*np.ones(phi_st.size)

# insert main field line start point (in MAG) as first entry 
phiMAG = MLON
thetaMAG = np.pi/2. - MLAT
R = 1.

u_st = np.insert(u_st, 0, R*np.sin(thetaMAG)*np.sin(phiMAG))
v_st = np.insert(v_st, 0, R*np.sin(thetaMAG)*np.cos(phiMAG))
w_st = np.insert(w_st, 0, R*np.cos(thetaMAG))

# Convert field line start points points from MAG to GSM
x_st = (np.nan)*np.empty((Nb+1,))
y_st = (np.nan)*np.empty((Nb+1,))
z_st = (np.nan)*np.empty((Nb+1,))
for i in range(Nb+1):
    v = ps.MAGtoGSM([u_st[i], v_st[i], w_st[i]], month, day, year, UT)
    x_st[i] = v[0]
    y_st[i] = v[1]
    z_st[i] = v[2]

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
for i in range(Nb+1):
    X0 = [x_st[i], y_st[i], z_st[i]] # Initial condition
    sol = odeint(dXds, X0, s_grid)
    solns[:,:,i] = sol

# initialize vectors for defining field line cut plane
v1=(np.nan)*np.empty((3,))
v2=(np.nan)*np.empty((3,))
v3=(np.nan)*np.empty((3,))
U1=(np.nan)*np.empty((3,))
U2=(np.nan)*np.empty((3,))
U3=(np.nan)*np.empty((3,))


# restrict the field lines to stop when reaching 1*R_E from the origin
solns_restr=[] # initialize list of np_arrays, one for each restricted field line
for i in range(Nb+1):  # loop over field lines
    # define condition on the field line points
    tr = np.logical_and(solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 >=1.,solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 < 20.)
    # create the arrays of the restricted field line componentwise
    solx=solns[:,0,i]
    solx=solx[tr]
    soly=solns[:,1,i]
    soly=soly[tr]
    solz=solns[:,2,i]
    solz=solz[tr]
    # reasemble and add to the list
    sol=np.column_stack([solx,soly,solz])
    solns_restr.append(sol)
    if (i == 0): 
        # do for main field line
        v1 = sol[0,:]
        v2 = sol[-1,:]
        v3 = sol[10,:]
        # define cut plane coordinates based on main field line 
        # (U3 normal to the plane)
        U2 = (v1-v2)/np.linalg.norm(v1-v2)
        U3 = np.cross(v3-v1, U2)/np.linalg.norm(np.cross(v3-v1, U2))
        U1 = np.cross(U2, U3)   

#    return [U1, U2, U3]

x_1d = np.linspace(0, 4, n)
y_1d = np.linspace(-3, 3, m)
X, Y = np.meshgrid(x_1d, y_1d) # grid of points on the cutplane
Z = np.zeros((n, m))
for i in range(n):
    for j in range(m):
        # grid of the corresponding values of variable. To be color plotted
        Z[i,j]=data_in_U(parameter,X[i,j],Y[i,j],U1,U2)

#------------------------------
kameleon.close()
print("Closed " + filename)
#-------------------------------

from_list=solns_restr[0]
sol=np.array(from_list)

f = open(fname,'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('A dataset with one polyline and no attributes\n')
f.write('ASCII\n')
f.write('\n')
f.write('DATASET POLYDATA\n')
f.write('POINTS '+str(sol.shape[0])+' float\n')
for k in range(sol.shape[0]):
    f.write('%e %e %e\n'%(sol[k,0],sol[k,1],sol[k,2]))

f.write('LINES '+'1'+' '+str(sol.shape[0]+1)+'\n' )
f.write(str(sol.shape[0])+'\n')
for k in range(sol.shape[0]):
    f.write(str(k)+'\n')

f.close()
