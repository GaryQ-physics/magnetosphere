# field_line_vtk_script

import sys
import os
import numpy as np
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')
sys.path.append(conf["m_path"] + 'magnetosphere/events/')
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
debug = False

# run parameters
Nlong = 5
Nb = 6
sign=-1  # changes sign of magnetic field used to trace the field lines

# Plot title
title = 'SCARR5 ' + str(year) + '-' + str(month) + '-' + str(day) + 'T07:00'
filename = conf["f_path"] + '3d__var_3_e' + str(year) + str(month) + str(day) + '-070000-000.out.cdf'

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

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
    if( x**2 + y**2 + z**2 >=1.):
        return data
    else:
        return 0.

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
    if 1e-9<Bm<1e+7:
        return (sign/Bm)*B
    else:
        if(debug):
            if(Bm >= 1e+7): print('FIELD TOO HIGH')
            if(Bm <= 1e-7): print('FIELD TOO LOW')
        return [0., 0., 0.] #or nan


u_st = (np.nan)*np.empty((Nb+1,))
v_st = (np.nan)*np.empty((Nb+1,))
w_st = (np.nan)*np.empty((Nb+1,))

R = 1.01
eps=3.*deg
for i in range(Nb+1):
    if(i==0):
        delta=0.
    elif(i>Nb/2):
        #delta=(i-Nb/2)*eps
        delta=(-i)*eps
    else:
        #delta=(i-Nb/2)*eps
        delta=(-i)*eps
    phi = MLON
    theta = np.pi/2. - MLAT + delta
    u_st[i] = R*np.sin(theta)*np.sin(phi)
    v_st[i] = R*np.sin(theta)*np.cos(phi)
    w_st[i] = R*np.cos(theta)
    

# Convert field line start points points from MAG to GSM
x_st = (np.nan)*np.empty((Nb+1,))
y_st = (np.nan)*np.empty((Nb+1,))
z_st = (np.nan)*np.empty((Nb+1,))
for i in range(Nb+1):
    v = ps.MAGtoGSM([u_st[i], v_st[i], w_st[i]], month, day, year, UT)
    x_st[i] = v[0]
    y_st[i] = v[1]
    z_st[i] = v[2]

# Trace field lines
s_grid = np.linspace(0, 200., 2000.)
if(debug): print(s_grid)
if(debug): print(s_grid.size)
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
    #tr = np.logical_and(solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 >=1.,solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 <= 224.**2+128.**2+128.**2)
    # create the arrays of the restricted field line componentwise
    went_out=0
    came_back=0
    end_val=solns.shape[0]-1
    for k in range(solns.shape[0]):
        if(solns[k,0,i]**2 + solns[k,1,i]**2 + solns[k,2,i]**2 > 1.1**2): went_out=1
        if(solns[k,0,i]**2 + solns[k,1,i]**2 + solns[k,2,i]**2 < 1.1**2 and went_out==1):
            came_back=1
            end_val=k
            break
    tr=np.arange(end_val+1)
    solx=solns[:,0,i]
    solx=solx[tr]
    soly=solns[:,1,i]
    soly=soly[tr]
    solz=solns[:,2,i]
    solz=solz[tr]
    # reasemble and add to the list
    sol=np.column_stack([solx,soly,solz])
    solns_restr.append(sol)
    if(debug and i==Nb+1): print(solns[:,:,i])
    if(debug and i==Nb+1): print(sol)

#------------------------------
kameleon.close()
print("Closed " + filename)
#-------------------------------

mlong_array=[0., 10.*deg, -10.*deg, 20.*deg, -20.*deg]
for i in range(Nb+1+Nlong):
    out_fname=conf["m_path"] + 'magnetosphere/data/' + 'field_line'+str(i)+'.vtk'
    if(i > Nb):
        mlong=mlong_array[i-Nb-1]
        mlat=np.linspace(-np.pi/2,np.pi/2,100)
        sol=np.column_stack([np.cos(mlong)*np.cos(mlat), np.sin(mlong)*np.cos(mlat), np.sin(mlat)])
        print('writing ' + out_fname)
        f = open(out_fname,'w')
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
        print('closed ' + out_fname)
    else:
        from_list=solns_restr[i]
        sol=np.array(from_list)
        print('writing ' + out_fname)
        f = open(out_fname,'w')
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
        print('closed ' + out_fname)
        f.close()
