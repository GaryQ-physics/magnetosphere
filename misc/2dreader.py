#from spacepy import pybats as pb
import spacepy.pybats.bats as bt
import numpy as np
from matplotlib import pyplot as plt
############################################################################

filename = "y=0_var_1_e20191012-020200-016.out"


# create the pybats object which contains our 2D magnetosphere data
data2d = bt.Bats2d(filename)

# to list the keys within data2d:
print(data2d.keys())

# I want to extract p, bx, by, bz
p, bx, by, bz =  np.array(data2d['p']), np.array(data2d['bx']), np.array(data2d['by']), np.array(data2d['bz'])

# now plot them:


plt.clf()

# plot p first:
ax1 = plt.subplot(121)
ax1.set_aspect('equal', 'box')
data2d.add_contour('x', 'z', 'p', target = ax1)

plt.grid(True)

ax2 = plt.subplot(122)
ax2.set_aspect('equal', 'box')

plt.xlim([-20, 20])
plt.ylim([-20, 20])

data2d.add_contour('x', 'z', 'bz', target = ax2, dolog = True, cmap = 'jet')

# add open and closed magnetic field lines:
#data2d.add_b_magsphere(target = ax2)

plt.grid(True)
plt.show()


