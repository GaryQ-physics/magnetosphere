# tester_for_ps

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()

sys.path.append(conf["m_path"] + 'magnetosphere/events/')
import pos_sun as ps

# [year, month, day, hours, minutes, seconds, x, y, z, 'car']
# or
# [year, month, day, hours, minutes, seconds, r, long in degrees, latitude in degrees, 'sph']

data = [ 
        [2003.0, 11.0, 20.0, 7.0,   0.0, 0.0, 0.,   0. ,  1.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 0.,   0.,   1.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,  68.5, 50.0, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 43.0, 0.0, 1., 166.0, 52.5, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 47.0, 0.0, 1., 166.0, 52.5, 'sph']
    ]

# https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi
# GSM [x, y, z, MLT]
sscweb = [
            [-0.46, -0.00, 0.89,  2.64457],
            [-0.17,  0.00, 0.99, 13.16680],
            [-0.09,  0.64, 0.76, 17.73347],
            [-0.72, -0.18, 0.67,  1.13969],
            [-0.72, -0.19, 0.67,  1.20329]
        ]

k = 0
for d in data:
    print('------------------------------------------------------------')
    syst = d[9] # sph or cart
    time = d[0:6]
    if syst == 'car':
        pos = d[6:9]
        v = ps.MAGtoGSM(pos, time, 'car', 'car')
        MLT = ps.MLTfromMAG(pos, time)
    elif syst == 'sph':
        r, MLON, MLAT = d[6:9]
        v = ps.MAGtoGSM([r, MLAT, MLON], time, 'sph', 'car')
        # Compute MLT given GSM using equation 93 in
        # https://arxiv.org/abs/1611.10321
        MLT = ps.MLTfromMAG(MLON, time)
    else:
        print('INVALID COORDINATE TYPE. Use "car" or "sph"')
    UT = time[3] + time[4]/60.+ time[5]/3600.
    spacepy = np.append(v, MLT)
    print('Input:')
    print('   [Year, Month, Day, Hour, Minute, Second]')
    print('   [{0:.0f}, {1:.0f}, {2:.0f}, {3:.0f}, {4:.0f}, {5:.0f}]'.format(time[0], time[1], time[2], time[3], time[4], time[5]))
    print('   MAG           x       y      z')
    print('                {0:.2f}  {1:.2f} {2:.2f}'.format(d[6], d[7], d[8]))
    print('Output:')
    print('   GSM           x      y    z    MLT')
    print('   SpacePy:    {0:.2f}  {1:.2f} {2:.2f} {3:.2f}'.format(spacepy[0], spacepy[1], spacepy[2], spacepy[3]))
    print('   SSCWeb:     {0:.2f}  {1:.2f} {2:.2f} {3:.2f}'.format(sscweb[k][0], sscweb[k][1], sscweb[k][2], sscweb[k][3]))
    k = k + 1


