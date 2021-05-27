import numpy as np
from config import conf
from magnetometers import GetMagnetometerLocation
import util

def get_mhd_values(run, time, obs_point, var, norcut=False):
    if isinstance(obs_point,str):
        obs_point_str = obs_point
        if obs_point == "origin":
            obs_point = np.zeros(3)
        else:
            obs_point = GetMagnetometerLocation(obs_point_str, time, 'GSM', 'car')
    else:
        obs_point_str = '[%f,%f,%f]'%(obs_point[0],obs_point[1],obs_point[2])

    namefromstitch = conf[run+'_derived']+'timeseries/slices/' \
                        + var+'_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
                        + '_obs_point=%s'%(obs_point_str)

    if norcut:
        namefromstitch = namefromstitch + '.npy'
    else:
        namefromstitch = namefromstitch  \
              + '_rcut=%f.npy'%(util.get_rCurrents(run))

    return np.load(namefromstitch)
