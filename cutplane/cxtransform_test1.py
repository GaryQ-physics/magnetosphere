import os
import sys
import numpy as np

# Add path of config.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import cxtransform as cx

# [year, month, day, hours, minutes, seconds, x, y, z, 'car']
# or
# [year, month, day, hours, minutes, seconds, r, long in degrees, latitude in degrees, 'sph']

data = [ 
        [2003.0, 11.0, 20.0, 7.0,   0.0, 0.0, 1.,   0. ,  0.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,   0.,   0.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,  68.5, 50.0, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 43.0, 0.0, 1., 166.0, 52.5, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 47.0, 0.0, 1., 166.0, 52.5, 'sph']
    ]
############ 
# before the first two were:
# data=[
#        [2003.0, 11.0, 20.0, 7.0,   0.0, 0.0, 0.,   0. ,  1.,  'car'],
#        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 0.,   0.,   1.,  'car'],
# with corresponding:
# sscweb=[
#          [-0.46, -0.00, 0.89,  2.64457],
#          [-0.17,  0.00, 0.99, 13.16680],

# This was a poor choice since cartesian 0, 0, 1 is at a singularity of polar coordinates so phi not well defined.
# coincidentally two things were true: 
#    1) There was a bug in cx.MAGtoMLT where when passed a (3,) array, it interpereted it as three phi values and returned array of three MLT's for each one
#       and only the the first one would appear in the print statement of this test. (the first being phi = 0 since 0, 0, 1 passed)
#    2) The sscweb apparently uses the convention that the zenith 0, 0, 1 has 0 longitude
# thus the MLT's agreed by accident 

# It just so happens np.arctan2(0., 0.) == 0 and sscweb was using the same convention aparently.
# So phi was set to zero anyway as it would be for cartesian 1,0,0, which is why the MLT's agree
######################


# Results from
# https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi
#   unfortunately can only manually enter spherical coordinates for MAG

# GSM [x, y, z, MLT]
sscweb = [
            [-0.68, -0.64, -0.36,  2.64457],
            [ 0.94,  0.30,  0.16, 13.16680],
            [-0.09,  0.64,  0.76, 17.73347],
            [-0.72, -0.18,  0.67,  1.13969],
            [-0.72, -0.19,  0.67,  1.20329]
        ]


k = 0
for d in data:
    print('------------------------------------------------------------')
    syst = d[9] # sph or cart
    time = d[0:6]
    if syst == 'car':
        pos = d[6:9]
        v = cx.MAGtoGSM(pos, time, 'car', 'car')
        MLT = cx.MAGtoMLT(pos, time)
    elif syst == 'sph':
        r, MLON, MLAT = d[6:9]
        v = cx.MAGtoGSM([r, MLAT, MLON], time, 'sph', 'car')
        MLT = cx.MAGtoMLT(MLON, time, onlyMLON = True)
    else:
        print('INVALID COORDINATE TYPE. Use "car" or "sph"')
    UT = time[3] + time[4]/60.+ time[5]/3600.

    spacepy = np.append(v, MLT) # NOTE: if MLT has more then one element, this will be longer than 4 and only 
                                # the first element of MLT will show up in the print statements that follow

    #print(spacepy.shape)

    print('Input:')
    print('   [Year, Month, Day, Hour, Minute, Second]')
    print('   [{0:.0f}, {1:.0f}, {2:.0f}, {3:.0f}, {4:.0f}, {5:.0f}]'.format(time[0], time[1], time[2], time[3], time[4], time[5]))
    if syst == 'car':
        print('   MAG           x       y      z')
    else:
        print('   MAG           r     mlon    mlat')
    print('                {0:.2f}  {1:.2f} {2:.2f}'.format(d[6], d[7], d[8]))
    print('Output:')
    print('   GSM           x      y    z    MLT')
    print('   SpacePy:    {0:.2f}  {1:.2f}  {2:.2f}  {3:.2f}'.format(spacepy[0], spacepy[1], spacepy[2], spacepy[3]))
    print('   SSCWeb:     {0:.2f}  {1:.2f}  {2:.2f}  {3:.2f}'.format(sscweb[k][0], sscweb[k][1], sscweb[k][2], sscweb[k][3]))
    print('   Diff:       {0:.3f} {1:.3f} {2:.3f} {3:.3f}'.format(sscweb[k][0] - spacepy[0],\
                                                                   sscweb[k][1] - spacepy[1],\
                                                                   sscweb[k][2] - spacepy[2],\
                                                                   sscweb[k][3] - spacepy[3]))
    k = k + 1
