# L_shell

import sys
sys.path.append('../../../../lib/python2.7/site-packages/')
import kameleon_pull as kp
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from scipy.integrate import odeint

filename='/home/gary/3d__var_3_e20031120-070000-000.out.cdf'

kp.k_open(filename)

def Bx(x,y,z):
	ret=kp.main([filename,'bx',x,y,z])
	return ret

def By(x,y,z):
	ret=kp.main([filename,'by',x,y,z])
	return ret
def Bz(x,y,z):
	ret=kp.main([filename,'bz',x,y,z])
	return ret

phi_st=np.linspace(0, 2*np.pi, 10)
theta=np.pi/8.
y_st=np.sin(theta)*np.sin(phi_st)
z_st=np.sin(theta)*np.cos(phi_st)
x_st=np.cos(theta)*np.ones(phi_st.size)

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
	return [Bx(U[0],U[1],U[2])/L, By(U[0],U[1],U[2])/L, Bz(U[0],U[1],U[2])/L]

t_even = np.linspace(0, 10., 100)
solns = (np.nan)*np.empty((t_even.size, 3, x_st.size))
for i in range(x_st.size):
	U0 = [x_st[i], y_st[i], z_st[i]]    # Initial condition
	sol = odeint(dUdt, U0, t_even)  # method='RK23', t_eval=None)
	solns[:,:,i]=sol

#------------------------------
kp.k_close()
#-------------------------------






fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax = plt.axes(projection='3d')

#for k in range(sol.shape[0]):
#    ax.scatter(sol[k, 0], sol[k, 1], sol[k, 0], marker=m)
for i in range(x_st.size):
	ax.plot3D(solns[:,0,i], solns[:,1,i], solns[:,2,i], 'gray')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.set_xlim(-3,3)
ax.set_ylim(-3,3)
ax.set_zlim(-3,3)
#plt.plot(sol[:, 0], sol[:, 1], '.')

plt.show()
