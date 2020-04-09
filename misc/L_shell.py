# L_shell

import sys
sys.path.append('../../../../lib/python2.7/site-packages/')
#import kameleon_pull as kp
sys.path.append('../../../../lib/python2.7/site-packages/ccmc/')
import _CCMC as ccmc
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pos_sun as ct
#coordinate_transformations_magnetosphere           pos_sun


filename='/home/gary/3d__var_3_e20031120-070000-000.out.cdf'
month=11
day=20
year=2003
kameleon = ccmc.Kameleon()
kameleon.open(filename)
print(filename, "OPENED")
interpolator = kameleon.createNewInterpolator()

def filePull(argv):
	#print("main ran")
	if (len(argv) == 4):
	    variable = argv[0]
	    c0 = float(argv[1])
	    c1 = float(argv[2])
	    c2 = float(argv[3])
	    kameleon.loadVariable(variable)
	    var = interpolator.interpolate(variable, c0, c1, c2)
	    return var
	else:
		print('Usage: <filename> <variable> x, y, z \n python kameleon_test rho -40 0 0')

def Bx(x,y,z):
	#ret=kp.main([filename,'bx',x,y,z])
	ret=filePull(['bx',x,y,z])
	return ret

def By(x,y,z):
	ret=filePull(['by',x,y,z])
	return ret
def Bz(x,y,z):
	ret=filePull(['bz',x,y,z])
	return ret

sign=+1

"""
phi_st=0.
theta=np.linspace(0, 0.5*np.pi, 10)
R=1.
u_st=R*np.sin(theta)*np.sin(phi_st)
v_st=R*np.sin(theta)*np.cos(phi_st)
w_st=R*np.cos(theta)

"""
phi_st=np.linspace(0, 2*np.pi, 30)
theta=np.pi/2.
R=3.
u_st=R*np.sin(theta)*np.sin(phi_st)
v_st=R*np.sin(theta)*np.cos(phi_st)
w_st=R*np.cos(theta)*np.ones(phi_st.size)


x_st=np.zeros((np.size(u_st),)) #initialize 
y_st=np.zeros((np.size(u_st),))
z_st=np.zeros((np.size(u_st),))
for i in range(np.size(u_st)):
	v=ct.MAGtoGSM([u_st[i],v_st[i],w_st[i]],month,day,year)
	print(v)
	x_st[i]=v[0]
	y_st[i]=v[1]
	z_st[i]=v[2]

"""
phi_st=np.linspace(0, 2*np.pi, 30)
theta=np.pi/2.
R=3.
x_st=R*np.sin(theta)*np.sin(phi_st)
y_st=R*np.sin(theta)*np.cos(phi_st)
z_st=R*np.cos(theta)*np.ones(phi_st.size)


phi_st=0.
theta=np.linspace(0, 0.5*np.pi, 10)
R=1.
x_st=R*np.sin(theta)*np.sin(phi_st)
y_st=R*np.sin(theta)*np.cos(phi_st)
z_st=R*np.cos(theta)
"""
print("---------")
print(u_st)
print(v_st)
print(w_st)
print("--------")
print(x_st)
print(y_st)
print(z_st)

X0=1.
Y0=0.
Z0=0.
n=50
m=50
x_1d = np.linspace(-3, 3, n)
y_1d = np.linspace(-3, 3, m)
x, y = np.meshgrid(x_1d, y_1d)

B_x = np.empty((n, m)) # Bx on  grid
B_x[:] = np.nan 
B_y = np.empty((n, m)) # By on grid
B_y[:] = np.nan 
B_z = np.empty((n, m)) # By on grid
B_z[:] = np.nan 
L=0
for i in range(n): # iterate over rows
    if(i==0): print("i=", i)
    if(i==1): print("i=", i)
    if(i==40): print("i=", i)
    for j in range(m): # iterate over columns    
        #x[i,j]=i*gridsize
        #y[i,j]=j*gridsize
        L=np.sqrt(Bx(x[i,j],y[i,j],0)*Bx(x[i,j],y[i,j],0)+By(x[i,j],y[i,j],0)*By(x[i,j],y[i,j],0))
        B_x[i,j]=Bx(x[i,j],y[i,j],0)*np.log(L)/L
        B_y[i,j]=By(x[i,j],y[i,j],0)*np.log(L)/L
        #print(i,j,x[i,j],y[i,j])

# dx/dt = Bx(x,y,z)
# dy/dt = By(x,y,z)
# dz/dt = By(x,y,z)
# Let U = [x, y]
# dU[0]/dt = Bx(U[0],U[1],U[2])
# dU[1]/dt = By(U[0],U[1],U[2])
# dU[2]/dt = Bz(U[0],U[1],U[2])

# New function
def dUdt(U, t):
	L=np.sqrt(Bx(U[0],U[1],U[2])*Bx(U[0],U[1],U[2])+By(U[0],U[1],U[2])*By(U[0],U[1],U[2])+Bz(U[0],U[1],U[2])*Bz(U[0],U[1],U[2]))
	if 1e-6<L<1e+6:
		return [sign*Bx(U[0],U[1],U[2])/L, sign*By(U[0],U[1],U[2])/L, sign*Bz(U[0],U[1],U[2])/L]
	else:
		return [0., 0., 0.] #or nan

t_even = np.linspace(0, 10., 100)
solns = (np.nan)*np.empty((t_even.size, 3, x_st.size))
for i in range(x_st.size):
	U0 = [x_st[i], y_st[i], z_st[i]]    # Initial condition
	sol = odeint(dUdt, U0, t_even)  # method='RK23', t_eval=None)
	solns[:,:,i]=sol

#------------------------------
#kp.k_close()
kameleon.close()
print("KAMELEON CLOSED")
#-------------------------------






fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax = plt.axes(projection='3d')

#for k in range(sol.shape[0]):
#    ax.scatter(sol[k, 0], sol[k, 1], sol[k, 0], marker=m)
for i in range(x_st.size):
	ax.plot3D(solns[:,0,i], solns[:,1,i], solns[:,2,i], 'gray')
ax.plot3D(x_st,y_st,z_st,'red')
ax.plot3D(u_st,v_st,w_st,'green')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.set_xlim(-5,5)
ax.set_ylim(-5,5)
ax.set_zlim(-5,5)
#plt.plot(sol[:, 0], sol[:, 1], '.')

plt.show()