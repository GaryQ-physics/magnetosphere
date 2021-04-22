import numpy as np
from numba import njit
import pandas as pd
import datetime
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
from swmf_file_reader.read_swmf_files import read_all, find_index, interpolate
from derivatives import get_partials
from magnetometers import GetMagnetometerLocation
from named_var_indexes import _x,_y,_z,_b1x,_b1y,_b1z
import pandas as pd

@njit
def _jit_probe_b1():
    # needed ?
    pass


def slice_probe_b1(run, time, obs_point, cache=None):
    ftag = util.time2CDFfilename(run,time)[:-8]
    if cache is None:
        cache = read_all(ftag)

    if isinstance(obs_point,str):
        obs_point_str = obs_point
        if obs_point == "origin":
            obs_point = np.zeros(3)
        else:
            obs_point = GetMagnetometerLocation(obs_point_str, time, 'GSM', 'car')
    else:
        obs_point_str = '[%f,%f,%f]'%(obs_point[0],obs_point[1],obs_point[2])

    outname = conf[run+'_derived'] + 'timeseries/slices/' \
        + 'b1_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
        + '_obs_point=%s.npy'%(obs_point_str)

    b1x = interpolate(ftag, obs_point, var='b1x', cache=cache)
    b1y = interpolate(ftag, obs_point, var='b1y', cache=cache)
    b1z = interpolate(ftag, obs_point, var='b1z', cache=cache)

    ########### check things for obs_point being a naive gridpoint
    indx = find_index(ftag, obs_point, cache=cache)
    if not obs_point_str in ['origin','colaba']:
        assert(indx is not None)

    if indx is not None:
        assert(obs_point[0] == cache['DataArray'][(0,*indx)])
        assert(obs_point[1] == cache['DataArray'][(1,*indx)])
        assert(obs_point[2] == cache['DataArray'][(2,*indx)])

        tx = np.allclose( b1x , cache['DataArray'][(_b1x,*indx)] )
        ty = np.allclose( b1y , cache['DataArray'][(_b1y,*indx)] )
        tz = np.allclose( b1z , cache['DataArray'][(_b1z,*indx)] )
        if not (tx and ty and tx):
            print(b1x , cache['DataArray'][(_b1x,*indx)])
            print(b1y , cache['DataArray'][(_b1y,*indx)])
            print(b1z , cache['DataArray'][(_b1z,*indx)])
    ###########

    np.save(outname, np.array([b1x,b1y,b1z], dtype=np.float32))


def stitch_probe_b1(run, times, obs_point):
    if isinstance(obs_point,str):
        obs_point_str = obs_point
    else:
        obs_point_str = '[%f,%f,%f]'%(obs_point[0],obs_point[1],obs_point[2])

    columns = ['b1x','b1y','b1z']
    df_name = conf[run+'_derived']+'timeseries/df_b1' \
        + '_obs_point=%s.pkl'%(obs_point_str)

    dtimes = []
    slice_arrays = []
    for time in list(times):
        outname = conf[run+'_derived'] + 'timeseries/slices/' \
            + 'b1_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
            + '_obs_point=%s.npy'%(obs_point_str)

        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        if os.path.exists(outname):
            slice_arrays.append(np.load(outname))
        else:
            raise ValueError ('no data for time '\
                +str(time)+'. \nrun slice_probe_b1(...)')

    df = pd.DataFrame(data=slice_arrays, columns = columns,
                        index=dtimes)
    df.to_pickle(df_name)

