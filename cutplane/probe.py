import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

from util import time2filename, filemeta
import _CCMC as ccmc

def probe(time, P, var=None, debug=False):

    if type(time) == str:
        filename = time
    else:
        filename = time2filename(time)

    if not os.path.exists(filename):
        raise ValueError('Not found: ' + filename)
        return

    meta = filemeta(filename)
    kameleon = ccmc.Kameleon()

    kameleon.open(filename)

    interpolator = kameleon.createNewInterpolator()

    ret = {}
    if var is None:
        # Get data for all parameters
        for key in meta['parameters']:
            kameleon.loadVariable(key)
            ret[key] = interpolator.interpolate(key, P[0], P[1], P[2])
    elif type(var) == str:
        kameleon.loadVariable(var)
        ret = interpolator.interpolate(var, P[0], P[1], P[2])
    else:
        if type(var) != list and type(var) != tuple:
            raise ValueError('var must be None, str, list, or tuple')
    
        for v in var:
            if type(v) != str:
                raise ValueError('var must be a str, or a list/tuple of strs')              

            kameleon.loadVariable(v)
            ret[v] = interpolator.interpolate(v, P[0], P[1], P[2])
            
    kameleon.close()

    return ret