import numpy as np
import cxtransform as cx


#magnetometer_stations = {'YKClat':62.480, 'YKClon':245.518,
#                         'FRNlat':37.0913, 'FRNlon':-119.7193,
#                         'FURlat':48.17, 'FURlon':11.28 }


magnetometer_GEOs = {'YKC' : (62.480, 245.518), # (lat, lon) in GEO
                'FRN' : (37.0913, -119.7193),
                'FUR' : (48.17, 11.28),
                'colaba' : (18.907, 72.815) }

def GetMagnetometerLocation(station, time, csys_out, ctype_out):
    lat, lon = magnetometer_GEOs[station]
    return cx.transform(np.array([1., lat, lon]), time, 'GEO', csys_out, ctype_in='sph', ctype_out=ctype_out)

