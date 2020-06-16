# Installation (did not work with python 3.7.7)
# conda create -n kamodo python==3.7 # Double equal is important!
# conda activate kamodo
# pip install kamodo
# pip install decorify
# pip install execnet
# Download a Kameleon-formatted CDF file, e.g.,
# http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070102-000.out.cdf
# Edit paths

# See also https://github.com/nasa/Kamodo
# Documentation for kameleon wrapper
# https://ccmc.gsfc.nasa.gov/Kamodo/notebooks/kameleon-kamodo/

if True:
    python_path = '/Users/robertweigel/kameleon/bin/python'
    fname = '/Users/robertweigel/git/magnetosphere/data/SCARR5_GM_IO2/IO2/'
    fname = fname + '3d__var_3_e20031120-070102-000.out.cdf'

if False:
    python_path = '/home/gary/magnetosphere/kameleon/bin/python'
    fname = '/home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/'
    fname = fname + '3d__var_3_e20031120-070000-000.out.cdf'

import time
import numpy as np
from kamodo.readers._kameleon_kamodo import Kameleon

kameleon = Kameleon(fname, python_path, 'p')
points = np.linspace(2,10,999999).reshape((-1,3))

to = time.time()
data = kameleon.p(points)
tf = time.time()
print('{0:d} points in {1:.1f} s = {2:.1f} pts/s' \
      .format(data.shape[0], tf-to, data.shape[0]/(tf-to))) 
# 19.5 s    
