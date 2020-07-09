import os
import sys
import numpy as np

# Add path of config.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import cxtransform as cx

data_car = np.array([ 
        [2003.0, 11.0, 20.0, 7.0,   0.0, 0.0, 1.,   0. ,  0.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,   0.,   0.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,   0.,   0.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,   0.,   0.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,   0.,   0.,  'car']
    ])

data_sph = np.array([ 
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,  68.5, 50.0, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 43.0, 0.0, 1., 166.0, 52.5, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 47.0, 0.0, 1., 166.0, 52.5, 'sph']
    ])


time = np.array(data_car[:,0:6], dtype=float)
pos = np.array(data_car[:,6:9], dtype=float)

v = cx.MAGtoGSM(pos, time, 'car', 'car')
MLT = cx.MAGtoMLT(pos, time, csys='car')

for i in range(data_car.shape[0]):
    print( str(v[i,:]) + ' : ' + str(cx.MAGtoGSM(pos[i,:], time[i,:], 'car', 'car')) )
    print( str(MLT[i]) + ' : ' + str(cx.MAGtoMLT(pos[i,:], time[i,:], csys='car')) )

time = np.array(data_sph[:,0:6], dtype=float)
pos = np.array(data_sph[:,6:9], dtype=float)

v = cx.MAGtoGSM(pos, time, 'sph', 'car')
MLT = cx.MAGtoMLT(pos[:,1], time)

for i in range(data_sph.shape[0]):
    print( str(v[i,:]) + ' : ' + str(cx.MAGtoGSM(pos[i,:], time[i,:], 'sph', 'car')) )
    print( str(MLT[i]) + ' : ' + str(cx.MAGtoMLT(pos[i,1], time[i,:])) )
