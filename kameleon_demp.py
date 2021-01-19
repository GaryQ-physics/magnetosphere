#script

import os
import sys
import numpy as np

execfile('config.py')
import _CCMC as ccmc

#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/3d__var_1_t00002000_n0002973.out.cdf'


#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/IONO-2D.swmf.it061215_130700_000.cdf' # 405.5 kB
#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/IONO-2D.swmf.it061215_060700_000.cdf' # 2.7 MB
filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf' # 2.7 MB

kameleon = ccmc.Kameleon()
print 0
kameleon.open(filename)
print 1

numvars = kameleon.getNumberOfVariables()
print(numvars)

print 2

interpolator = kameleon.createNewInterpolator()
kameleon.loadVariable('x')
kameleon.loadVariable('y')
kameleon.loadVariable('z')

out = interpolator.interpolate('y', 10., 10., 10.)
print(out)

print 3

A = np.linspace(0,100,10)
B = np.linspace(0,100,10)

'''

for i in range(A.size):
    #print((i,j))
    for n in range(1000):
        #print(n)
        a = interpolator.interpolate('x', 10.*np.random.rand(), A[i], np.random.rand())
        b = interpolator.interpolate('x', 10.*np.random.rand(), A[i], np.random.rand())
        print(a)
        print(b)
        print(np.abs(a-b))
        a = interpolator.interpolate('y', 10.*np.random.rand(), A[i], np.random.rand())
        b = interpolator.interpolate('y', 10.*np.random.rand(), A[i], np.random.rand())
        assert(a==b)
        a = interpolator.interpolate('z', 10.*np.random.rand(), A[i], np.random.rand())
        b = interpolator.interpolate('z', 10.*np.random.rand(), A[i], np.random.rand())
        assert(a==b)
        #assert(np.abs(a-b) < 1e-9)
'''

for i in range(A.size):
    for j in range(B.size):
        #print((i,j))
        for n in range(1000):
            #print(n)
            a = interpolator.interpolate('x', 10.*np.random.rand(), A[i], B[j])
            b = interpolator.interpolate('x', 10.*np.random.rand(), A[i], B[j])
            assert(a==b)
            a = interpolator.interpolate('y', 10.*np.random.rand(), A[i], B[j])
            b = interpolator.interpolate('y', 10.*np.random.rand(), A[i], B[j])
            assert(a==b)
            a = interpolator.interpolate('z', 10.*np.random.rand(), A[i], B[j])
            b = interpolator.interpolate('z', 10.*np.random.rand(), A[i], B[j])
            assert(a==b)
            #assert(np.abs(a-b) < 1e-9)

print 4



import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

print(sys.version)
print(sys.path)
print(os.environ['PYTHONPATH'].split(os.pathsep))


fig = plt.figure()
# https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html
ax = plt.axes(projection='3d') 

sol0 = np.array([[0,0,0],[1,1,1],[2,2,2],[3,3,3]])
sol1 = np.array([[0,0,0],[-1,1,1],[-2,2,2],[-3,3,3]])

print(sol0.shape)
print(sol1.shape)

#ax.plot3D(sol0[:,0], sol0[:,1], sol0[:,2], 'blue', lw=2, alpha=0.7)

for i in range(1):

#for i in range(A.size):
    sol = np.zeros((B.size,3))
    for j in range(B.size):
        sol[j,0] = interpolator.interpolate('x', np.random.rand(), A[i], B[j])
        sol[j,1] = interpolator.interpolate('y', np.random.rand(), A[i], B[j])
        sol[j,2] = interpolator.interpolate('z', np.random.rand(), A[i], B[j])
    print(sol)
    ax.plot3D(sol[:,0], sol[:,1], sol[:,2], 'b', lw=2, alpha=0.7)


#ax.plot3D(sol1[:,0], sol1[:,1], sol1[:,2], 'blue', lw=2, alpha=0.7)


ax.set(xlabel = "X")
ax.set(ylabel = "Y")
ax.set(zlabel = "Z")

L=4
ax.set_xlim(-L, L)
ax.set_ylim(-L, L)
ax.set_zlim(-L, L)

plt.show()

#import matplotlib.pyplot as plt
#plt.plot(A, a)
#plt.show()
