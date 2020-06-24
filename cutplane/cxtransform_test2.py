import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import _CCMC as ccmc
import cxtransform as cx

# Inputs for file and magnetic field model
year = 2003
day = 20
month = 11
hours = 7
minutes = 0
seconds = 0
time = (year, month, day, hours, minutes, seconds)

filename = conf['run_path'] + '3d__var_3_e' \
                + str(year) + str(month) + str(day) + '-070000-000.out.cdf'

from datetime import datetime
epoch = datetime(2000, 1, 1, 12, 0, 0)
delta = datetime(year, month, day, hours, minutes, seconds) - epoch
EphemTime = int(delta.total_seconds())

# open kameleon
print("Opening " + filename)
kameleon = ccmc.Kameleon()
kameleon.open(filename)
print("Opened " + filename)
interpolator = kameleon.createNewInterpolator()

# no arguments assumes native coordinate system
coordinate_interpolator = kameleon.createCoordinateInterpolator() 
print('File epoch time: ' + str(coordinate_interpolator.getEphemTime()) + ' seconds')
coordinate_interpolator.setEphemTime(EphemTime)
print('Computed epoch time: ' + str(coordinate_interpolator.getEphemTime()) \
      + ' seconds since 2000-01-01T12:00:00')

n = 4 # Number of tests
MAG = 3.*np.random.rand(n, 3) # Random points to test

print('Model coordinates: ' + coordinate_interpolator.get_model_coords())
coordinate_interpolator.setPreferredCoordinates('MAG')

for i in range(n):
    MAG_p = ccmc.Position()
    MAG_p.c0 = MAG[i, 0]
    MAG_p.c1 = MAG[i, 1]
    MAG_p.c2 = MAG[i, 2]
    GSM_p = ccmc.Position()
    coordinate_interpolator.convertCoordinates(MAG_p, GSM_p)
    GSM_kameleon = [GSM_p.c0, GSM_p.c1, GSM_p.c2]

    GSM_spacepy = cx.MAGtoGSM([MAG[i,0], MAG[i,1], MAG[i,2]], time, 'car', 'car')
    
    print('------------------------------------------------')
    print('MAG           x        y        z')
    print('         {0:8.2f} {1:8.2f} {2:8.2f}'.format(MAG[i,0], MAG[i,1], MAG[i,2]))
    print('GSM           x        y        z')
    print('Kameleon {0:8.2f} {1:8.2f} {2:8.2f}'
          .format(GSM_kameleon[0], GSM_kameleon[1], GSM_kameleon[2]))
    print('SpacePy  {0:8.2f} {1:8.2f} {2:8.2f}'
          .format(GSM_spacepy[0], GSM_spacepy[1], GSM_spacepy[2]))
    print('Diff     {0:8.2e} {1:8.2e} {2:8.2e}'
          .format(GSM_spacepy[0] - GSM_kameleon[0],
                  GSM_spacepy[1] - GSM_kameleon[1],
                  GSM_spacepy[2] - GSM_kameleon[2]))

print('------------------------------------------------')
kameleon.close()
