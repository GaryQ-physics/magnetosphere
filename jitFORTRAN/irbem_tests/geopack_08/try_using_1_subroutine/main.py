import numpy as np
import pymod

try:
    print(pymod.geomag_08_v)
    print('yes pymod.geomag_08_v')
except:
    print('no pymod.geomag_08_v')

try:
    print(pymod.geomag_08)
    print('yes pymod.geomag_08')
except:
    print('no pymod.geomag_08')


N = 5

IYEAR = 1
IDAY = 1
IHOUR = 1
MIN = 1
ISEC = 1

J = 1

XGEO_V = np.array( [0., 1., 2., 3., 4.] , dtype=np.float64)
YGEO_V = np.array( [0., 1., 2., 3., 4.] , dtype=np.float64)
ZGEO_V = np.array( [0., 1., 2., 3., 4.] , dtype=np.float64)
XMAG_V = np.array( [0., 1., 2., 3., 4.] , dtype=np.float64)
YMAG_V = np.array( [0., 1., 2., 3., 4.] , dtype=np.float64)
ZMAG_V = np.array( [0., 1., 2., 3., 4.] , dtype=np.float64)

print(XGEO_V)
print(YGEO_V)
print(ZGEO_V)
print(XMAG_V)
print(YMAG_V)
print(ZMAG_V)

pymod.geomag_08_v(XGEO_V,YGEO_V,ZGEO_V,XMAG_V,YMAG_V,ZMAG_V,
                    J,
                    IYEAR,IDAY,IHOUR,MIN,ISEC,
                    N)

print(XGEO_V)
print(YGEO_V)
print(ZGEO_V)
print(XMAG_V)
print(YMAG_V)
print(ZMAG_V)

