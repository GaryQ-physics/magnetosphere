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

data=[ [2003.0, 11.0, 20.0, 7.0, 0.0, 0.0, 0., 0. ,1.,'car'],
    [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 0.0, 0.,1.,'car'],
    [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1., 68.5, 50.0,'sph'],
    [2003.0, 11.0, 20.0, 18.0, 43.0, 0.0, 1., 166.0, 52.5, 'sph'],
    [2003.0, 11.0, 20.0, 18.0, 47.0, 0.0, 1., 166.0, 52.5, 'sph'] ]

'''
# MAG inputs
# corresponding NASA GSM  : MLT

[2003.0, 11.0, 20.0, 7.0, 0.0, 0.0, 0., 0. ,1.,'car']
-0.46    -0.00    0.89    : 2.64457

[2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 0., 0.,1.,'car']
-0.17    0.00    0.99    : 13.16680

[2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1., 68.5, 50.0,'sph']
-0.09    0.64    0.76    : 17.73347

[2003.0, 11.0, 20.0, 18.0, 43.0, 0.0, 1., 166.0, 52.5, 'sph']
-0.72   -0.18    0.67    : 1.13969

[2003.0, 11.0, 20.0, 18.0, 47.0, 0.0, 1., 166.0, 52.5, 'sph']
-0.72   -0.19    0.67    : 1.20329
'''

for l in data:
    print '--------'
    year,month,day,hours,minutes,seconds,a,b,c,syst = l
    Time = [year,month,day,hours,minutes,seconds]
    if(syst=='car'):
        x,y,z = a,b,c
        v = ps.MAGtoGSM([x,y,z],Time,'car','car')
    elif(syst=='sph'):
        r,MLON,MLAT = a,b,c
        #x,y,z = ps.StoC( , (90.-MLAT)*deg, MLON*deg )
        v = ps.MAGtoGSM([r,MLAT,MLON],Time,'sph','car')
        phi=MLON*deg
        subsol_pt = ps.GSMtoMAG([1,0,0],Time,'car','car')
        phi0 = np.arctan2(subsol_pt[1], subsol_pt[0])
        delta = phi-phi0
        if(delta > np.pi):
            delta=delta-2.*np.pi
        elif(delta <= -np.pi):
            delta=delta+2.*np.pi
        MLT = 12.*hr + delta*24.*hr/(360*deg)
        print MLT
    else:
        print 'INVALID COORDINATE TYPE'
    UT=hours*hr + minutes*minn + seconds*s
    #print v[3]==hours
    #print v[4]==minutes
    #print v[5]==seconds
    #print v[0],v[1],v[2]
    print v
