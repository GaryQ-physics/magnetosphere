# plane_cut


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
def p(x,y,z):
	ret=filePull(['p',x,y,z])
	return ret

deg=(np.pi/180)

sign=-1

phi_st=np.linspace(0, 2*np.pi, 30)
theta=np.pi/2.
R=3.
u_st=R*np.sin(theta)*np.sin(phi_st)
v_st=R*np.sin(theta)*np.cos(phi_st)
w_st=R*np.cos(theta)*np.ones(phi_st.size)

MLON=68.50*deg  #68.50
MLAT=50.00*deg  #50.00
phiMAG=MLON
thetaMAG=np.pi/2. - MLAT
R=1.
u_st=np.insert(u_st,0,R*np.sin(thetaMAG)*np.sin(phiMAG))
v_st=np.insert(v_st,0,R*np.sin(thetaMAG)*np.cos(phiMAG))
w_st=np.insert(w_st,0,R*np.cos(thetaMAG))

x_st=np.zeros((np.size(u_st),)) #initialize 
y_st=np.zeros((np.size(u_st),))
z_st=np.zeros((np.size(u_st),))
for i in range(np.size(u_st)):
	v=ct.MAGtoGSM([u_st[i],v_st[i],w_st[i]],month,day,year)
	print(v)
	x_st[i]=v[0]
	y_st[i]=v[1]
	z_st[i]=v[2]



print("---------")
print(u_st)
print(v_st)
print(w_st)
print("--------")
print(x_st)
print(y_st)
print(z_st)


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


plt.clf()
fig = plt.figure(figsize=plt.figaspect(0.5))
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax2 = fig.add_subplot(1,2,2)

v1=(np.nan)*np.empty((3,))
v2=(np.nan)*np.empty((3,))
v3=(np.nan)*np.empty((3,))
U1=(np.nan)*np.empty((3,))
U2=(np.nan)*np.empty((3,))
U3=(np.nan)*np.empty((3,))
for i in range(u_st.size):
	tr=np.array( solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 >=1. )
	solx=solns[:,0,i]
	solx=solx[tr]
	soly=solns[:,1,i]
	soly=soly[tr]
	solz=solns[:,2,i]
	solz=solz[tr]
	sol=np.column_stack([solx,soly,solz])
	if (i==0):
		ax1.plot3D(sol[:,0], sol[:,1], sol[:,2], 'red')
		v1=sol[0,:]
		v2=sol[-1,:]
		v3=sol[20,:]
		U2 = (v1-v2)/np.linalg.norm(v1-v2)
		U3 = np.cross(v3-v1,U2)/np.linalg.norm(np.cross(v3-v1,U2))
		U1 = np.cross(U2,U3)
		solCut=np.zeros((sol.shape[0],2))
		for k in range(sol.shape[0]):
			solCut[k,0]=np.dot(sol[k,:],U1)
			solCut[k,1]=np.dot(sol[k,:],U2)
		print(solCut)
		ax2.plot(solCut[:,0],solCut[:,1], 'red')

	else:
		ax1.plot3D(sol[:,0], sol[:,1], sol[:,2], 'gray')

def BU1(u,v):
	x,y,z =u*U1+v*U2
	B=np.array([Bx(x,y,z),By(x,y,z),Bz(x,y,z)])
	return np.dot(B,U1)
def BU2(u,v):
	x,y,z =u*U1+v*U2
	B=np.array([Bx(x,y,z),By(x,y,z),Bz(x,y,z)])
	return np.dot(B,U2)
def BU3(u,v):
	x,y,z =u*U1+v*U2
	B=np.array([Bx(x,y,z),By(x,y,z),Bz(x,y,z)])
	return np.dot(B,U3)
def p_U(u,v):
	x,y,z =u*U1+v*U2
	return p(x,y,z)

n=50
m=50
x_1d = np.linspace(0, 4, n)
y_1d = np.linspace(-3, 3, m)
X, Y = np.meshgrid(x_1d, y_1d)
Z=np.zeros((n,m))
for i in range(n):
	for j in range(m):
		Z[i,j]=p_U(X[i,j],Y[i,j])

ax2.pcolormesh(X, Y, Z)



#------------------------------
#kp.k_close()
kameleon.close()
print("KAMELEON CLOSED")
#-------------------------------

para1, para2 = np.meshgrid(np.linspace(-1., 3., n), np.linspace(-2., 2., m))
X_slice=U1[0]*para1+U2[0]*para2+v1[0]*np.ones((n,m))
Y_slice=U1[1]*para1+U2[1]*para2+v1[1]*np.ones((n,m))
Z_slice=U1[2]*para1+U2[2]*para2+v1[2]*np.ones((n,m))

X_o=-1*para1+0*para2+v1[0]*np.ones((n,m))
Y_o=0*para1+0*para2+v1[1]*np.ones((n,m))
Z_o=0*para1+1*para2+v1[2]*np.ones((n,m))

ax1.plot_surface(X_slice, Y_slice, Z_slice)
ax1.plot_surface(X_o, Y_o, Z_o)

#ax.set(xlabel='x', ylabel='y', zlabel='z')


plt.show()
