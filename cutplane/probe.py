import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import numpy as np
from util import time2filename, filemeta
import _CCMC as ccmc

def probe(time, P, var=None, debug=False, dictionary=False):

    P = np.array(P)
    if P.shape == (3, ):
        P = np.array([P])

    if type(time) == str:
        filename = time
    else:
        filename = time2filename(time)

    if not os.path.exists(filename):
        raise ValueError('Not found: ' + filename)

    kameleon = ccmc.Kameleon()
    kameleon.open(filename)
    interpolator = kameleon.createNewInterpolator()

    if var is None:
        meta = filemeta(filename)
        # Get data for all parameters, store in dictionary
        ret = {}
        for key in meta['parameters']:
            kameleon.loadVariable(key)
            arr = np.nan*np.empty((P.shape[0],))
            for k in range(P.shape[0]):
                arr[k] = interpolator.interpolate(key, P[k,0], P[k,1], P[k,2])
            if arr.size == 1:
                arr = arr[0]
            ret[key] = arr

    elif type(var) == str:
        kameleon.loadVariable(var)
        ret = np.nan*np.empty((P.shape[0],))
        for k in range(P.shape[0]):
            ret[k] = interpolator.interpolate(var, P[k,0], P[k,1], P[k,2])
        if ret.size == 1:
            ret = ret[0]

    else:
        if type(var) != list and type(var) != tuple:
            raise ValueError('var must be None, str, list, or tuple')
        if dictionary:
            ret = {}
        else:
            ret = np.nan*np.empty((P.shape[0], len(var)))
        i = 0
        for v in var:
            if type(v) != str:
                raise ValueError('var must be a str, or a list/tuple of strs')              

            kameleon.loadVariable(v)
            arr = np.nan*np.empty((P.shape[0],))
            for k in range(P.shape[0]):
                arr[k] = interpolator.interpolate(v, P[k,0], P[k,1], P[k,2])
            if arr.size == 1:
                arr = arr[0]
            if dictionary:
                ret[v] = arr
            else:
                ret[:,i] = arr
                i = i + 1
        if not dictionary:
            if P.shape[0] == 1:
                ret = ret.flatten()

    kameleon.close()

    return ret
