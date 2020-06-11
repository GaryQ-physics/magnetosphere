import sys
import os
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import _CCMC as ccmc

#fname = '/Users/robertweigel/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070102-000.out.cdf'
fname = '/home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070000-000.out.cdf'


points = np.linspace(2,10,9999999).reshape((-1,3))

print(points)
print(points.shape)

# open kameleon 
kameleon = ccmc.Kameleon()
kameleon.open(fname)
print("Opened " + fname)
interpolator = kameleon.createNewInterpolator()

a = np.nan*np.empty((3333333,))
for i in range(3333333):
    kameleon.loadVariable('rho')
    a[i] = interpolator.interpolate('rho', points[i,0], points[i,1], points[i,2])

# close kameleon 
kameleon.close()
print("Closed " + fname)

print(a)
print(a.shape)
