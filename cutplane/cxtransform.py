import sys
import os
import numpy as np

# Add path of config.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import spacepy.coordinates as sc
from spacepy.time import Ticktock




def MAGtoGSM(v_MAG, time, ctype_in, ctype_out):
    
    in_type = type(v_MAG)    
    v_MAG = np.array(v_MAG)
    cvals = sc.Coords(v_MAG, 'MAG', ctype_in)
    cvals.ticks = Ticktock(tstr(time), 'ISO')
    newcoord = cvals.convert('GSM', ctype_out)
    v_GSM = newcoord.data[0,:]
    
    if in_type == np.ndarray:
        return v_GSM
    else:
        return v_GSM.tolist()


def GSMtoMAG(v_GSM, time, ctype_in, ctype_out):
    """Convert from GSM to MAG coordinates using SpacePy library"""

    in_type = type(v_GSM)
    v_GSM = np.array(v_GSM)
    cvals = sc.Coords(v_GSM, 'GSM', ctype_in)
    cvals.ticks = Ticktock(tstr(time), 'ISO')
    newcoord = cvals.convert('MAG', ctype_out)
    v_MAG = newcoord.data[0,:]
    
    if in_type == np.ndarray:
        return v_MAG
    else:
        return v_MAG.tolist()

def tstr(time):
    """Create ISO 8601 date/time string given array of integers
    
    tstr((2000, 1, 1)) # 2000-01-01T00:00:00
    tstr((2000, 1, 1, 2)) # 2000-01-01T02:00:00
    tstr((2000, 1, 1, 2, 3)) # 2000-01-01T02:03:00
    tstr((2000, 1, 1, 2, 3, 4)) # 2000-01-01T02:03:04
    """
    
    # TODO: Use datetime function for formatting
    
    time = np.array(time)
    if len(time.shape) > 1:
        ret = []
        for i in range(time.shape[0]):
            ret.append(tstr(time[i]))
        return ret
    
    time = list(time)
    # TODO: Check that time is valid
    if len(time) < 4:
        # TODO: Throw error
        pass
    
    if len(time) > 6:
        time = time[0:6]
    else:
        pad = 6 - len(time)
        time = time + pad*[0]
    
    return '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(time)


def transform(v, time, csys_in, csys_out, ctype_in=None, ctype_out=None):
    """Transfrom between coordinates systems using SpacePy library.

    Parameters
    ----------
    v : array-like
        list of three floats
        list containing lists of three floats
        np.array of three floats
        (Nv, 3) float np.array, where Nv=1 or Nv=Nt
        
    time : array-like
        list of 3+ ints
        list containing lists of 3+ ints
        np.array of 3+ ints
        (Nt, 3) float np.array, where Nt = 1 or Nt = Nv

        The 3+ ints are [year, month, day, [hours, [minutes, [seconds]]]]

    csys_in : str

    csys_out : str
    
    ctype_in : str
               'car' (default) or 'sph'

    ctype_out : str
               'car' (default) or 'sph'
               
    Returns
    -------
    array-like with structure matching either `time` (if `Nt` != 1) or
    `v` (if `Nv` =! 1). If `Nv` and `Nt` != 1, structure matches `v`. 

    Examples
    --------
    t1 = [2000, 1, 1, 0, 0, 0]
    t2 = [2000, 1, 1, 2, 0, 0]
    v1 = [0, 0, 1]
    v2 = [0, 1, 0]

    # All equivalent and return a list with three floats    
    transform([0, 1, 1], time1)
    transform([0, 1, 1], time1, ctype_in='car')
    transform([0, 1, 1], time1, ctype_in='car', ctype_out='car')

    # The following 3 calls return a list with two lists of 3 elements
    # 1. Transform two vectors at same time t1
    transform([v1, v2], t1) 

    # 2. Transform two vectors at different times
    transform([v1, v2], [t1, t2])

    # 3. Transform one vector at different times
    transform(v1, [t1, t2])

    # The following returns a (2, 3) np.array
    # Transform one vector at three times
    transform(np.array([v1, v2]), t1) 

    """
    
    import numpy.matlib
    
    if ctype_in is None:
        ctype_in = 'car'
    if ctype_out is None:
        ctype_out = 'car'

    in_type = type(v)
    v = np.array(v)
    t = np.array(time)

    if len(t.shape) > 1 and len(v.shape) > 1:
        if t.shape[1] != v.shape[1]:
            raise ValueError("")
    if len(t.shape) > 1:
        v = numpy.matlib.repmat(v, t.shape[0], 1)
    if len(t.shape) == 1:
        t = numpy.matlib.repmat(t, v.shape[0], 1)

    cvals = sc.Coords(v, csys_in, ctype_in)
    
    cvals.ticks = Ticktock(tstr(t), 'ISO')
    newcoord = cvals.convert(csys_out, ctype_out)

    if in_type == np.ndarray:
        return newcoord.data
    else:
        return newcoord.data.tolist()

    
def GEOtoGSM(v_GEO, time, ctype_in, ctype_out):
    """Transform from GEO to GSM coordinates using SpacePy library
    
    """

    # return transform(v, time, 'GEO', 'GSM', ctype_in, ctype_out)
    
    in_type = type(v_GEO)

    v_GEO = np.array(v_GEO)
    cvals = sc.Coords(v_GEO, 'GEO', ctype_in)
    
    if v_GEO.shape == (3, ):
        cvals.ticks = Ticktock(tstr(time), 'ISO')
        newcoord = cvals.convert('GSM', ctype_out)
        v_GSM = newcoord.data[0, :]        
    else:
        # Create as many ticks as there are points in v_GEO
        cvals.ticks = Ticktock([tstr(time)]*v_GEO.shape[0], 'ISO')
        newcoord = cvals.convert('GSM', ctype_out)
        v_GSM = newcoord.data

    if in_type == np.ndarray:
        return v_GSM
    else:
        return v_GSM.tolist()


def StoC(r, theta, phi):
    """Convert from spherical to cartesian coordinates

    r, theta, phi: array-like

    """
    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)
    return x, y, z


def UTtoHMS(UT, **kwargs):
    """Convert universal time in fractional hours into integer hour, minutes, seconds
    
    from cxtransform import UTtoHMS

    print(UTtoHMS(12))              # [12, 0, 0]

    print(UTtoHMS(24))              # [0, 0, 0]

    print(UTtoHMS(24, keep24=True)) # [24, 0, 0]
    """

    keep24 = False
    if 'keep24' in kwargs:
        keep24 = kwargs['keep24']
        
    if UT > 24 or UT < 0:
        raise ValueError('Required: 0 <= UT <= 24.')

    hours = int(UT)
    minutes = int((UT-hours)*60.)
    seconds = int(round((UT-hours-minutes/60.)*3600.))
    if seconds == 60:
        seconds = 0
        minutes = minutes + 1
    if minutes == 60:
        minutes = 0
        hours = hours + 1

    if hours == 24 and keep24 == False:
        return [0, 0, 0]

    return [hours, minutes, seconds]


def MAGtoMLT(pos, time):
    """Compute magnetic local time given a UT and MAG position or longitude

    Uses equation 93 in https://arxiv.org/abs/1611.10321

    Usage: 
    ------
    import cxtransform as cx
    mlt = cx.MAGtoMLT(MAGlong, time)
    mlt = cx.MAGtoMLT([MAGx, MAGy, MAGz], time)

    Returns:
    --------
    mlt: float

    Example:
    --------
    import cxtransform as cx
    mlt = cx.MAGtoMLT(0., [2000, 1, 1, 0, 0, 0])
    print(mlt)

"""

    if isinstance(pos, float):
        phi = pos*np.pi/180
        pos = [pos]
    else:
        phi = np.arctan2(pos[1], pos[0])

    pos = np.array(pos)
    time = np.array(time)
    if isinstance(pos[0], np.ndarray):
        phi = np.arctan2(pos[:, 1], pos[:, 0])
    else:
        phi = pos*np.pi/180

    subsol_pt = transform(np.array([1, 0, 0]), time, 'GSM', 'MAG')
    #import pdb;pdb.set_trace()
    phi_cds = np.arctan2(subsol_pt[:, 1], subsol_pt[:, 0])
        
    delta = phi - phi_cds
    
    idx = np.where(delta > np.pi)
    delta[idx] = delta[idx] - 2.*np.pi

    idx = np.where(delta <= -np.pi)
    delta[idx] = delta[idx] + 2.*np.pi

    MLT = 12. + delta*24./(2.*np.pi)
    return MLT
