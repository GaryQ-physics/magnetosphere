import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/')
sys.path.append(conf["k_path"] + 'kameleon/lib/python2.7/site-packages/ccmc/')
import _CCMC as ccmc

year = 2003
day = 20
month = 11
hours = 7.
minutes = 0.
seconds = 0.

filename = conf["f_path"] + '3d__var_3_e' + str(year) + str(month) + str(day) + '-070000-000.out.cdf'

kameleon = ccmc.Kameleon()
kameleon.open(filename)
interpolator = kameleon.createNewInterpolator()
coordinate_interpolator = kameleon.createCoordinateInterpolator() # no arguments assumes native
print 'epoch time:', coordinate_interpolator.getEphemTime(), 'seconds'
coordinate_interpolator.setEphemTime(0)

# Set date to 2000 001 00:00:00
# Set lat, lon, Re to 1, 1, 1 at
# https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi
# Click "GM"
coordinate_interpolator.setPreferredCoordinates("MAG")
point = 1.,0.,0. # MAG cartesian [Re]

query_point = ccmc.Position()
query_point.c0, query_point.c1, query_point.c2 = point

model_point = ccmc.Position()
coordinate_interpolator.convertCoordinates(query_point,model_point)
print coordinate_interpolator.get_model_coords(),
print model_point.c0, model_point.c1, model_point.c2 # GSM
