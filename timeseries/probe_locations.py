import numpy as np
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
from swmf_file_reader.named_var_indexes import str2index,_b1x,_b1y,_b1z
import pandas as pd


def slice_probe(run, time, obs_point, var, cache=None):
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
    obs_point = np.array(obs_point)

    outname = conf[run+'_derived'] + 'timeseries/slices/' \
        + var+'_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
        + '_obs_point=%s.npy'%(obs_point_str)

    if var in ['b1','b','j','u']:
        Fx = interpolate(ftag, obs_point, var=var+'x', cache=cache)
        Fy = interpolate(ftag, obs_point, var=var+'y', cache=cache)
        Fz = interpolate(ftag, obs_point, var=var+'z', cache=cache)

        ########### check things for obs_point being a naive gridpoint
        if var == 'b1':
            indx = find_index(ftag, obs_point, cache=cache)
            if not obs_point_str in ['origin','colaba']:
                assert(indx is not None)

            if indx is not None:
                assert(obs_point[0] == cache['DataArray'][(0,*indx)])
                assert(obs_point[1] == cache['DataArray'][(1,*indx)])
                assert(obs_point[2] == cache['DataArray'][(2,*indx)])

                tx = np.allclose( Fx , cache['DataArray'][(_b1x,*indx)] )
                ty = np.allclose( Fy , cache['DataArray'][(_b1y,*indx)] )
                tz = np.allclose( Fz , cache['DataArray'][(_b1z,*indx)] )
                if not (tx and ty and tx):
                    print('WARNING')
                    print(Fx , cache['DataArray'][(_b1x,*indx)])
                    print(Fy , cache['DataArray'][(_b1y,*indx)])
                    print(Fz , cache['DataArray'][(_b1z,*indx)])
        ###########

        np.save(outname, np.array([Fx,Fy,Fz], dtype=np.float32))
    else:
        F = interpolate(ftag, obs_point, var=var, cache=cache)
        np.save(outname, np.array([F], dtype=np.float32))


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

