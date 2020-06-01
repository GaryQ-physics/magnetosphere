import sys
import os
import numpy as np

# Add path of config_paths.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import spacepy.coordinates as sc
from spacepy.time import Ticktock


def MAGtoGSM(v_MAG, time, ctype_in, ctype_out):
    """Convert from MAG to GSM coordinates using SpacePy library.

    MAG: array-like
    
    time: list or tuple of ints
          [year,month,day,hours,minutes,seconds]
    
    ctype_in: str
        'car' or 'sph'

    ctype_out: str
        'car' or 'sph'
          
    """
    
    v_MAG = np.array(v_MAG)
    cvals = sc.Coords(v_MAG, 'MAG', ctype_in)
    T = tuple(time)
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % T 
    cvals.ticks = Ticktock(t_str, 'ISO')
    newcoord = cvals.convert('GSM', ctype_out)
    v_GSM = newcoord.data[0,:]
    return v_GSM


def GSMtoMAG(v_GSM, time, ctype_in, ctype_out):
    """Convert from GSM to MAG coordinates using SpacePy library"""

    v_GSM = np.array(v_GSM)
    cvals = sc.Coords(v_GSM, 'GSM', ctype_in)
    T = tuple(time)
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % T
    cvals.ticks = Ticktock(t_str, 'ISO') # add ticks
    newcoord = cvals.convert('MAG', ctype_out)
    v_MAG = newcoord.data[0,:]
    return v_MAG

def GEOtoGSM(v_GEO, time, ctype_in, ctype_out):
    """Convert from GEO to GSM coordinates using SpacePy library
    
    Example:
    --------
    time = [2000, 1, 1, 0, 0, 0]
    GEOtoGSM([0, 1, 1], time, 'car', 'car')
    GEOtoGSM([[0, 1, 1],[0, 1, 1]], time, 'car', 'car')
    GEOtoGSM(np.array([[0, 1, 1],[0, 1, 1]]), time, 'car', 'car')
    
    """

    T = tuple(time)
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % T

    v_GEO = np.array(v_GEO)
    cvals = sc.Coords(v_GEO, 'GEO', ctype_in)
    
    if v_GEO.shape == (3, ):
        cvals.ticks = Ticktock(t_str, 'ISO') # add ticks
        newcoord = cvals.convert('GSM', ctype_out)
        v_GSM = newcoord.data[0,:]
        return v_GSM
    else:
        cvals.ticks = Ticktock([t_str]*v_GEO.shape[0], 'ISO') # add ticks
        newcoord = cvals.convert('GSM', ctype_out)
        v_GSM = newcoord.data
        return v_GSM


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
    
    from pos_sun import UTtoHMS
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
    """Compute magnetic local time given a universial time given MAG position or longitude

    Uses equation 93 in https://arxiv.org/abs/1611.10321

    Usage: 
    ------
    import pos_sun as ps
    mlt = ps.MAGtoMLT(MAGlong, time)
    mlt = ps.MAGtoMLT([MAGx, MAGy, MAGz], time)

    Returns:
    --------
    mlt: float

    Example:
    --------
    mlt = MAGtoMLT(0., [2000, 1, 0, 0, 0, 0])
    print(mlt)

"""

    if isinstance(pos, float):
        phi = pos*np.pi/180
    else:
        phi = np.arctan2(pos[1], pos[0])
    subsol_pt = GSMtoMAG([1, 0, 0], time, 'car', 'car')
    phi_cds = np.arctan2(subsol_pt[1], subsol_pt[0])
    delta = phi - phi_cds
    if delta > np.pi:
        delta = delta - 2.*np.pi
    elif delta <= -np.pi:
        delta = delta + 2.*np.pi
    MLT = 12. + delta*24./(2.*np.pi)
    return MLT
