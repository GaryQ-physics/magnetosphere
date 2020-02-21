# L_shell

import sys
sys.path.append('../../../../lib/python2.7/site-packages/')
import kameleon_pull as kp
import numpy as np
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

X0=1
Y0=0
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

U0 = [X0, Y0, Z0]    # Initial condition
t_even = np.linspace(0, 1., 100)
sol = odeint(dUdt, U0, t_even)  # method='RK23', t_eval=None)

#------------------------------
kp.k_close()
#-------------------------------

plt.plot(sol[:, 0], sol[:, 1], '.')

plt.show()
