# B_field_lines_write_demo1

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import B_field_lines_write as Bfl

print(sys.version)
print(sys.path)
print(os.environ['PYTHONPATH'].split(os.pathsep))

fig = plt.figure()
ax = plt.axes(projection='3d') #https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html

solns_original = Bfl.Compute([2003, 11, 20 , 7, 0, 0, 0, 176.00, 57.50], 6, use_grid=False)
solns_withGrid = Bfl.Compute([2003, 11, 20 , 7, 0, 0, 0, 176.00, 57.50], 6, use_grid=True)

print(len(solns_withGrid) == len(solns_original))

sol0 = solns_original[0]
sol1 = solns_withGrid[0]

for i in range(len(solns_withGrid)):
    sol0 = solns_original[i]
    sol1 = solns_withGrid[i]

    min_len = min(sol0.shape[0], sol1.shape[0])
    print(np.abs(sol1[:min_len,:] - sol0[:min_len,:]).max())

    print(sol0.shape)
    print(sol1.shape)
    print(sol0[0,:])
    print(sol1[0,:])

    ax.plot3D(sol0[:,0], sol0[:,1], sol0[:,2], 'red', lw=2)
    ax.plot3D(sol1[:,0], sol1[:,1], sol1[:,2], 'blue', lw=4, alpha=0.6)

ax.set(xlabel = "$X/R_E$ (GSM)")
ax.set(ylabel = "$Y/R_E$ (GSM)")
ax.set(zlabel = "$Z/R_E$ (GSM)")

L=20
ax.set_xlim(-L, L)
ax.set_ylim(-L, L)
ax.set_zlim(-L, L)

plt.show()
