# Run script with no arguments to see test results

import sys
import os
import numpy as np

# Add path of config_paths.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()

import spacepy.coordinates as sc
from spacepy.time import Ticktock

runTester = False

# units
deg = np.pi/180
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.


def MAGtoGSM(v_MAG, Time, ctype_in, ctype_out): # Time = [year,month,day,hours,minutes,seconds]
    '''Convert from MAG to GSM coordinates using SpacePy library'''
    
    if ( (ctype_in != 'sph' and ctype_in != 'car') or (ctype_out != 'sph' and ctype_out != 'car') ):
        print 'ctype ERROR'
        return (np.nan)*np.empty((3, ))
    v_MAG = np.array(v_MAG)
    cvals = sc.Coords(v_MAG, 'MAG', ctype_in)
    T = tuple(Time)
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % T # (year,month,day,hours,minutes,seconds)
    cvals.ticks = Ticktock(t_str, 'ISO') # add ticks
    newcoord = cvals.convert('GSM', ctype_out)
    v_GSM = np.array([newcoord.x[0],newcoord.y[0],newcoord.z[0]])
    if(runTester): v_GSM = np.array([newcoord.x[0],newcoord.y[0],newcoord.z[0],hours,minutes,seconds])
    return(v_GSM)

#function to convert coordinate vector in MAG coordinates to one in GSM coordinates
def GSMtoMAG(v_GSM,Time,ctype_in,ctype_out):
    if( (ctype_in != 'sph' and ctype_in != 'car') or (ctype_out != 'sph' and ctype_out != 'car') ):
        print 'ctype ERROR'
        return (np.nan)*np.empty((3,))
    v_GSM = np.array(v_GSM)
    cvals = sc.Coords(v_GSM, 'GSM', ctype_in)
    T=tuple(Time)
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % T # (year,month,day,hours,minutes,seconds)
    cvals.ticks = Ticktock(t_str, 'ISO') # add ticks
    newcoord = cvals.convert('MAG', ctype_out)
    v_MAG = np.array([newcoord.x[0],newcoord.y[0],newcoord.z[0]])
    if(runTester): v_MAG = np.array([newcoord.x[0],newcoord.y[0],newcoord.z[0],hours,minutes,seconds])
    return(v_MAG)

def GEOtoGSM(v_GEO,Time,ctype_in,ctype_out):
    if( (ctype_in != 'sph' and ctype_in != 'car') or (ctype_out != 'sph' and ctype_out != 'car') ):
        print 'ctype ERROR'
        return (np.nan)*np.empty((3,))
    v_GEO = np.array(v_GEO)
    cvals = sc.Coords(v_GEO, 'GEO', ctype_in)
    T=tuple(Time)
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % T # (year,month,day,hours,minutes,seconds)
    cvals.ticks = Ticktock(t_str, 'ISO') # add ticks
    newcoord = cvals.convert('GSM', ctype_out)
    v_GSM = np.array([newcoord.x[0],newcoord.y[0],newcoord.z[0]])
    return(v_GSM)

def GEOtoGSM_list(v_GEO,Time,ctype_in,ctype_out):
    if( (ctype_in != 'sph' and ctype_in != 'car') or (ctype_out != 'sph' and ctype_out != 'car') ):
        print 'ctype ERROR'
        return (np.nan)*np.empty((3,))
    cvals = sc.Coords(v_GEO, 'GEO', ctype_in)
    T=tuple(Time)
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % T # (year,month,day,hours,minutes,seconds)
    cvals.ticks = Ticktock([t_str]*len(v_GEO), 'ISO') # add ticks
    newcoord = cvals.convert('GSM', ctype_out)
    v_GSM = np.column_stack([newcoord.x,newcoord.y,newcoord.z])
    return(v_GSM)

#convert spherical to cartesian
def StoC(r,theta,phi):
    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)
    return x,y,z

def UTtoHMS(UT):
    hours = int(UT/hr)
    minutes = int((UT-hours*hr)/minn)
    seconds = int(round((UT-hours*hr-minutes*minn)/s))
    if(seconds==60):
        seconds=0
        minutes=minutes+1
    if(minutes==60):
        minutes=0
        hours=hours+1
    if(hours>=24):
        print 'WARNING: MIGHT BE NEXT DAY'
    return [hours,minutes,seconds]


if __name__ == "__main__":
    MAGtoGSM(v_MAG, Time, ctype_in, ctype_out)

