# Must be run in Python 2.7 (Kameleon Python wrapper requirement)
# Install Kameleon - see https://ccmc.gsfc.nasa.gov/Kameleon/Quick_start.html
# Download a Kameleon-formatted CDF file, e.g.,
# http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070102-000.out.cdf
# Edit paths

if True:
    ccmc_path = '/Users/robertweigel/kameleon/lib/python2.7/site-packages/ccmc/'
    fname = '/Users/robertweigel/git/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070102-000.out.cdf'

if False:
    ccmc_path = '/home/gary/magnetosphere/kameleon/lib/python2.7/site-packages/ccmc/'
    fname = '/home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070000-000.out.cdf'

import sys
import time
import numpy as np
sys.path.append(ccmc_path)
import _CCMC as ccmc

kameleon = ccmc.Kameleon()

kameleon.open(fname)
interpolator = kameleon.createNewInterpolator()
kameleon.loadVariable("p")
points = np.linspace(2,10,999999).reshape((-1,3))

i = 0
data = np.empty((points.shape[0]))
to = time.time()
for i in range(points.shape[0]):
    data[i] = interpolator.interpolate("p", points[i, 0], points[i, 1], points[i, 2])
kameleon.close()
tf = time.time()
print('{0:d} points in {1:.1f} s = {2:.1f} pts/s' \
      .format(data.shape[0], tf-to, data.shape[0]/(tf-to)))
#  1.6 s
