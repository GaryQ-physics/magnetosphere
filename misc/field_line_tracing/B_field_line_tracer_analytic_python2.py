# B_field_line_tracer_analytic

import sys
sys.path.append('../../../../lib/python2.7/site-packages/')
import numpy as np
import matplotlib.pyplot as plt
#import scipy as sp
from scipy.integrate import odeint
#import solve_ivp

def Bx(x,y,z):
	M = -31000.
	R=np.sqrt(x**2+y**2+z**2)
	ret=M*(3*x**2 - R**2)/R**5
	return ret
def By(x,y,z):
	M = -31000.
	R=np.sqrt(x**2+y**2+z**2)
	ret=3*M*y*x/R**5
	return ret
def Bz(x,y,z):
	M = -31000.
	R=np.sqrt(x**2+y**2+z**2)
	ret=3*M*z*x/R**5
	return ret

eps=0.0001
#gridsize=1./50.
#Delta_x=1.
#Delta_y=1.
l=10000
X0=0.2
Y0=0.1
n=50
m=50
#n=int(Delta_x/gridsize)
#m=int(Delta_y/gridsize)
x_1d = np.linspace(-1, 1, n)
y_1d = np.linspace(-1, 1, m)
x, y = np.meshgrid(x_1d, y_1d)
B_x = np.empty((n, m)) # Bx on  grid
B_x[:] = np.nan 
B_y = np.empty((n, m)) # By on grid
B_y[:] = np.nan 
s = np.empty((l,)) # parameter of line
s[:] = np.nan
X = np.empty((l,)) # x coord along line
X[:] = np.nan 
Y = np.empty((l,)) # y coord along line
Y[:] = np.nan 
i0=int(X0/eps)
j0=int(Y0/eps)
L=0
for i in range(n): # iterate over rows
    for j in range(m): # iterate over columns    
        #x[i,j]=i*gridsize
        #y[i,j]=j*gridsize
        L=np.sqrt(Bx(x[i,j],y[i,j],0)*Bx(x[i,j],y[i,j],0)+By(x[i,j],y[i,j],0)*By(x[i,j],y[i,j],0))
        B_x[i,j]=Bx(x[i,j],y[i,j],0)*np.log(L)/L
        B_y[i,j]=By(x[i,j],y[i,j],0)*np.log(L)/L

ds=eps
s[0]=0
X[0]=X0
Y[0]=Y0

for k in range(s.size-1):
	L=np.sqrt(Bx(X[k],Y[k],0)*Bx(X[k],Y[k],0)+By(X[k],Y[k],0)*By(X[k],Y[k],0))
	if L>0.000001:
		ds=eps/L
	else:
		ds=0
		print("B near ZERO")
	s[k+1] = s[k]+ds #determine parameterd
	X[k+1] = X[k]+ds*Bx(X[k],Y[k],0) #iterate eom
	Y[k+1] = Y[k]+ds*By(X[k],Y[k],0) #iterate eom
	#print(X[k],Y[k],s[k])


# dx/dt = Bx(x,y,0)
# dy/dt = By(x,y,0)
#
# Let U = [x, y]
# dU[0]/dt = Bx(U[0],U[1],0)
# dU[1]/dt = By(U[0],U[1],0)
 
# New function
def dUdt(U, t):
	Ux, Uy = U
	a=Bx(Ux,Uy,0.)
	b=By(Ux,Uy,0.)
	ret=[a, b]
	return ret
# New initial condition
U0 = [X0, Y0]    # Initial condition

t_even = np.linspace(0, 1.1e-8, 100)

"""
r = ode(dUdt, jac=None).set_integrator('zvode', method='bdf')
r.set_initial_value(U0, 0)
dt = 1000./30.
T=0
while r.successful() and r.t < 1000:
	T=T+dt
	print(r.t+dt, r.integrate(r.t+dt))
"""

sol = odeint(dUdt, U0, t_even)  #, method='RK23', t_eval=None)
print(sol[:,0])


fig, ax = plt.subplots()

q = ax.quiver(x, y, B_x, B_y, units='xy', angles='xy', scale=300)

plt.plot(X, Y)
plt.plot(sol[:, 0], sol[:, 1], '.')

plt.show()
