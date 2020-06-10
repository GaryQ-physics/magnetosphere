# Dependencies:
# pip install kamodo
# pip install execnet

from kamodo.readers._kameleon_kamodo import Kameleon
import numpy as np

python_path = '/Users/robertweigel/kameleon/bin/python'
fname = '/Users/robertweigel/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070102-000.out.cdf'
kameleon = Kameleon(fname, python_path, 'rho', 'p', 'bx','by','bz')

points = np.linspace(2,10,9999).reshape((-1,3))

# Interpolate rho onto points.
a = kameleon.rho(points)

print(a.shape)
