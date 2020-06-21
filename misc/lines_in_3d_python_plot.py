import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

print(sys.version)
print(sys.path)
print(os.environ['PYTHONPATH'].split(os.pathsep))


fig = plt.figure()
ax = plt.axes(projection='3d') #https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html

sol0 = np.array([[0,0,0],[1,1,1],[2,2,2],[3,3,3]])
sol1 = np.array([[0,0,0],[-1,1,1],[-2,2,2],[-3,3,3]])

print(sol0.shape)
print(sol1.shape)

ax.plot3D(sol0[:,0], sol0[:,1], sol0[:,2], 'red', lw=4)
ax.plot3D(sol1[:,0], sol1[:,1], sol1[:,2], 'blue', lw=4)

plt.show()
