# Dependencies:
# Python 3.7
# pip install kamodo
# pip install execnet

from kamodo.readers._kameleon_kamodo import Kameleon
import numpy as np

python_path = '/Users/robertweigel/kameleon/bin/python'
fname = '/Users/robertweigel/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070102-000.out.cdf'

kameleon = Kameleon(fname, python_path, 'p')

points = np.linspace(2,10,999999).reshape((-1,3))

print(points.shape)

# Interpolate rho onto points.
import time
to = time.time()
data = kameleon.p(points)
tf = time.time()
print('{0:d} points in {1:.1f} s = {2:.1f} pts/s' \
      .format(data.shape[0], tf-to, data.shape[0]/(tf-to))) 
    