# tester_for_ps

import sys
import os
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
sys.path.append(conf["m_path"] + 'magnetosphere/events/')
import pos_sun as ps

# units
deg=np.pi/180
amin=deg/60.
hr=1.
minn=hr/60.
s=minn/60.

data=[[2003.0, 11.0, 20.0, 7.0, 0.0, 0., 0. ,1.,'car'], [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 0.,1.,'car'],[2003.0, 11.0, 20.0, 17.0, 46.0, 1., 68.5, 50.0,'sph']]

'''
# MAG inputs
# corresponding NASA GSM  : MLT

[2003.0, 11.0, 20.0, 7.0, 0.0, 0., 0. ,1.,'car']
-0.46   -0.00    0.89   : 2.64457

[2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 0.,1.,'car']
-0.17    0.00    0.99    : 13.16680

[2003.0, 11.0, 20.0, 17.0, 46.0, 1., 68.5, 50.0,'sph']
-0.09    0.64    0.76   : 17.73347
'''

for l in data:
    seconds=0.
    year,month,day,hours,minutes,a,b,c,syst = l
    if(syst=='car'):
        x,y,z = a,b,c
    elif(syst=='sph'):
        r,MLON,MLAT = a,b,c
        x,y,z = ps.StoC( a, (90.-MLAT)*deg, MLON*deg )
    else:
        print 'INVALID SYST'
    #print day, z
    UT=hours*hr + minutes*minn + seconds*s
    v = ps.MAGtoGSM(np.array([x,y,z]),month,day,year,UT)
    print '--------'
    print v[3]==hours
    print v[4]==minutes
    print v[5]==seconds
    print v[0],v[1],v[2]
