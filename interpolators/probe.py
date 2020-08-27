import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

#from util import time2filename, filemeta
import util

'''
def interpolate(filename, variable, Q, library)
    if library == 'kameleonV':
        return kameleonV.interpolate(filename, Q[:,0], Q[:,1], Q[:,2], variable)

    if library == 'kameleon':
        kameleon = ccmc.Kameleon()
        kameleon.open(filename)
        interpolator = kameleon.createNewInterpolator()

        kameleon.loadVariable(variable)
        arr = np.nan*np.empty((P.variable[0],))
        for k in range(Q.shape[0]):
            arr[k] = interpolator.interpolate(variable, Q[k,0], Q[k,1], Q[k,2])
        return arr

    if library == 'pycdf'
        ii.interpolate(Q[:,0], Q[:,1], variable, filename)
'''

def probe(filename, P, var=None, debug=False, dictionary=False, library='kameleonV'):
    """
    library = 'kameleonV', 'kameleon', 'pycdf'
    """
    P = np.array(P)
    if P.shape == (3, ):
        P = np.array([P])

    assert(filename[0] == '/')
    #if type(time) == str:
    #    filename = filename
    #else:
    #    filename = util.time2filename(filename) #!!!!!!

    if not os.path.exists(filename):
        raise ValueError('Not found: ' + filename)

    ################### import apropriate file for library
    if library == 'kameleonV':
        assert(P.shape[1] == 3)
        sys.path.append(conf['interpolator'] + 'kameleonV-0.2.3/')
        import kameleonV
    if library == 'kameleon':
        assert(P.shape[1] == 3)
        sys.path.append(conf['interpolator'] + 'kameleon/lib/python2.7/site-packages/ccmc/')
        import _CCMC as ccmc
    if library == 'pycdf':
        assert(P.shape[1] == 2)
        sys.path.append(conf['interpolator'] + 'pycdf_with_scipy')
        import ionosphere_interpolator as ii
    ###################


    ################### define interpolate(variable, Q) funtion within this probe
    if library == 'kameleonV':
        def interpolate(variable, Q):
            return kameleonV.interpolate(filename, Q[:,0], Q[:,1], Q[:,2], variable)

    if library == 'kameleon':
        kameleon = ccmc.Kameleon()
        kameleon.open(filename)
        interpolator = kameleon.createNewInterpolator()

        def interpolate(variable, Q):
            kameleon.loadVariable(variable)
            arr = np.nan*np.empty((Q.shape[0],))
            for k in range(Q.shape[0]):
                arr[k] = interpolator.interpolate(variable, Q[k,0], Q[k,1], Q[k,2])
            return arr

    if library == 'pycdf':
        def interpolate(variable, Q):
            return ii.interpolate(Q[:,0], Q[:,1], variable, filename)
    ###################

    if var is None:
        meta = util.filemeta(filename)
        # Get data for all parameters, store in dictionary
        ret = {}
        for key in meta['parameters']:
            ret[key] = interpolate(key, P)

    elif type(var) == str:
        ret = interpolate(var, P)
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
            arr = interpolate(v, P)
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

    if library == 'kameleon':
        kameleon.close()

    return ret
