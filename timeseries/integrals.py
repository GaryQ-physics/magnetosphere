#python -c "import integrals as i; i.slice_B_biotsavart('DIPTSUR2',(2019,9,2,4,10,0,0),'colaba')"
#python -c "import integrals as i; i.stitch_B_biotsavart('DIPTSUR2',[(2019,9,2,4,10,0,0)],'colaba')"
import os
import numpy as np
from numba import njit
import pandas as pd
import datetime
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
from read_swmf_files2 import read_all
from derivatives import get_partials
from magnetometers import GetMagnetometerLocation
from named_var_indexes import _x,_y,_z

@njit
def _jit_B_biotsavart(DataArray, obs_point, rcut):
    unique_epsilons = np.array([  0.0625,
                                  0.1250,
                                  0.2500,
                                  0.5000,
                                  1.0000,
                                  2.0000,
                                  4.0000,
                                  8.0000 ])
    n_eps = unique_epsilons.size

    nBlock, nI, nJ, nK = DataArray.shape[1:]

    ret = np.zeros((n_eps,3), dtype=np.float32)
    for iBlockP in range(nBlock):
        epsilon = DataArray[_x,iBlockP,1,0,0] - DataArray[_x,iBlockP,0,0,0]
        assert(np.argwhere(unique_epsilons==epsilon).shape == (1,1))
        i_eps = np.argwhere(unique_epsilons==epsilon)[0][0]
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):
                    if i == 0 or j == 0 or k == 0 or i == nI-1 or j == nJ-1 or k == nK-1:
                        continue
                    distanceSquared = ( DataArray[_x, iBlockP, i,j,k]**2 \
                                      + DataArray[_y, iBlockP, i,j,k]**2 \
                                      + DataArray[_z, iBlockP, i,j,k]**2 )
                    if distanceSquared < rcut**2:
                        continue

                    integrand = np.empty((3,),dtype=np.float32); integrand[:]=np.nan

                    partials = get_partials(DataArray, 'b1', iBlockP,i,j,k)

                    curl_B1_tens = partials - partials.transpose()
                    curl_B1_x = curl_B1_tens[1,2]
                    curl_B1_y = curl_B1_tens[2,0]
                    curl_B1_z = curl_B1_tens[0,1]

                    r_x = obs_point[0] - DataArray[_x, iBlockP, i, j, k]
                    r_y = obs_point[1] - DataArray[_y, iBlockP, i, j, k]
                    r_z = obs_point[2] - DataArray[_z, iBlockP, i, j, k]
                    r = np.sqrt(r_x**2 + r_y**2 + r_z**2)
                    if r < 1e-5:
                        continue

                    integrand[2] = curl_B1_x*r_y - curl_B1_y*r_x
                    integrand[0] = curl_B1_y*r_z - curl_B1_z*r_y
                    integrand[1] = curl_B1_z*r_x - curl_B1_x*r_z
                    integrand[:] =  integrand[:]/r**3

                    ret[i_eps,:] = ret[i_eps,:] + integrand

    for i_eps in range(n_eps):
        ret[i_eps,:] = ( (unique_epsilons[i_eps]**3)/(4*np.pi) ) * ret[i_eps,:]
    return ret


def slice_B_biotsavart(run, time, obs_point, rcut=None, cache=None):
    '''
    inputs:
      run:
        string. run name
      time:
        1d np array. A single time of the simulation run
      obs_point:
        (3,) array like or string. The observation point at which the biot savart
        integral is evaluated for. If array, then its the xyz GSM coords. Otherwise,
        its a keyword like "colaba" or "origin". Note, something like "colaba"
        may have different gsm coordinates at different times
      rcut=None:
        float. The inner radius boundary, all points less the rcut from
        the origin will not contribute to the integral.
      cache=None:
        optional dictionary of cached data for that time slice. If None,
        then it is generated by reading the swmf data files.
    '''
    if rcut is None:
        rcut = util.get_rCurrents(run)

    if cache is None:
        cache = read_all(util.time2CDFfilename(run,time)[:-8])

    if isinstance(obs_point,str):
        obs_point_str = obs_point
        if obs_point == "origin":
            obs_point = np.zeros(3)
        else:
            obs_point = GetMagnetometerLocation(obs_point_str, time, 'GSM', 'car')
    else:
        obs_point_str = '[%f,%f,%f]'%(obs_point[0],obs_point[1],obs_point[2])

    outname = conf[run+'_derived'] + 'timeseries/slices/' \
        + 'B_biotsavart_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
        + '_obs_point=%s_rcut=%f.npy'%(obs_point_str, rcut)

    ret = _jit_B_biotsavart(cache['DataArray'], np.array(obs_point), rcut)
    np.save(outname+'-byepsilons.npy', ret)
    np.save(outname, np.sum(ret, axis=0))


def stitch_B_biotsavart(run, times, obs_point, rcut=None):
    '''
    inputs:
      run:
        string. run name
      times:
        2d array like of integers. Each row being the time of a timeslice.
      obs_point:
        (3,) numpy array. The observation point at which the biot savart
        integral is evaluated for.
      rcut=None:
        float. The inner radius boundary, all points less the rcut from
        the origin will not contribute to the integral.
    '''
    if rcut is None:
        rcut = util.get_rCurrents(run)

    if isinstance(obs_point,str):
        obs_point_str = obs_point
    else:
        obs_point_str = '[%f,%f,%f]'%(obs_point[0],obs_point[1],obs_point[2])

    columns = ['B_biotsavart_x','B_biotsavart_y','B_biotsavart_z']
    df_name = conf[run+'_derived']+'timeseries/df_biotsavart' \
        + '_obs_point=%s_rcut=%f.pkl'%(obs_point_str,rcut)

    dtimes = []
    slice_arrays = []
    for time in list(times):
        outname = conf[run+'_derived'] + 'timeseries/slices/' \
            + 'B_biotsavart_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
            + '_obs_point=%s_rcut=%f.npy'%(obs_point_str,rcut)

        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        if os.path.exists(outname):
            slice_arrays.append(np.load(outname))
        else:
            raise ValueError ('no data for time '\
                +str(time)+'. \nrun slice_B_biotsavart(...)')

    df = pd.DataFrame(data=slice_arrays, columns = columns,
                        index=dtimes)
    df.to_pickle(df_name)


@njit
def _jit_B_coulomb(DataArray, obs_point, rcut):
    unique_epsilons = np.array([  0.0625,
                                  0.1250,
                                  0.2500,
                                  0.5000,
                                  1.0000,
                                  2.0000,
                                  4.0000,
                                  8.0000 ])
    n_eps = unique_epsilons.size

    nBlock, nI, nJ, nK = DataArray.shape[1:]

    ret = np.zeros((n_eps,3), dtype=np.float32)
    for iBlockP in range(nBlock):
        epsilon = DataArray[_x,iBlockP,1,0,0] - DataArray[_x,iBlockP,0,0,0]
        assert(np.argwhere(unique_epsilons==epsilon).shape == (1,1))
        i_eps = np.argwhere(unique_epsilons==epsilon)[0][0]
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):
                    if i == 0 or j == 0 or k == 0 or i == nI-1 or j == nJ-1 or k == nK-1:
                        continue
                    distanceSquared = ( DataArray[_x, iBlockP, i,j,k]**2 \
                                      + DataArray[_y, iBlockP, i,j,k]**2 \
                                      + DataArray[_z, iBlockP, i,j,k]**2 )
                    if distanceSquared < rcut**2:
                        continue

                    integrand = np.empty((3,),dtype=np.float32); integrand[:]=np.nan

                    partials = get_partials(DataArray, 'b1', iBlockP,i,j,k)

                    div_B1 = partials[0,0] + partials[1,1] + partials[2,2]

                    r_x = obs_point[0] - DataArray[_x, iBlockP, i, j, k]
                    r_y = obs_point[1] - DataArray[_y, iBlockP, i, j, k]
                    r_z = obs_point[2] - DataArray[_z, iBlockP, i, j, k]
                    r = np.sqrt(r_x**2 + r_y**2 + r_z**2)
                    if r < 1e-5:
                        continue

                    integrand[0] = div_B1*r_x
                    integrand[1] = div_B1*r_y
                    integrand[2] = div_B1*r_z
                    integrand[:] =  integrand[:]/r**3

                    ret[i_eps,:] = ret[i_eps,:] + integrand

    for i_eps in range(n_eps):
        ret[i_eps,:] = ( (unique_epsilons[i_eps]**3)/(4*np.pi) ) * ret[i_eps,:]
    return ret


def slice_B_coulomb(run, time, obs_point, rcut=None, cache=None):
    '''
    inputs:
      run:
        string. run name
      time:
        1d np array. A single time of the simulation run
      obs_point:
        (3,) array like or string. The observation point at which the  coulomb
        integral is evaluated for. If array, then its the xyz GSM coords. Otherwise,
        its a keyword like "colaba" or "origin". Note, something like "colaba"
        may have different gsm coordinates at different times
      rcut=None:
        float. The inner radius boundary, all points less the rcut from
        the origin will not contribute to the integral.
      cache=None:
        optional dictionary of cached data for that time slice. If None,
        then it is generated by reading the swmf data files.
    '''
    if rcut is None:
        rcut = util.get_rCurrents(run)

    if cache is None:
        cache = read_all(util.time2CDFfilename(run,time)[:-8])

    if isinstance(obs_point,str):
        obs_point_str = obs_point
        if obs_point == "origin":
            obs_point = np.zeros(3)
        else:
            obs_point = GetMagnetometerLocation(obs_point_str, time, 'GSM', 'car')
    else:
        obs_point_str = '[%f,%f,%f]'%(obs_point[0],obs_point[1],obs_point[2])

    outname = conf[run+'_derived'] + 'timeseries/slices/' \
        + 'B_coulomb_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
        + '_obs_point=%s_rcut=%f.npy'%(obs_point_str, rcut)

    ret = _jit_B_coulomb(cache['DataArray'], np.array(obs_point), rcut)
    np.save(outname+'-byepsilons.npy', ret)
    np.save(outname, np.sum(ret, axis=0))


def stitch_B_coulomb(run, times, obs_point, rcut=None):
    '''
    inputs:
      run:
        string. run name
      times:
        2d array like of integers. Each row being the time of a timeslice.
      obs_point:
        (3,) numpy array. The observation point at which the coulomb
        integral is evaluated for.
      rcut=None:
        float. The inner radius boundary, all points less the rcut from
        the origin will not contribute to the integral.
    '''
    if rcut is None:
        rcut = util.get_rCurrents(run)

    if isinstance(obs_point,str):
        obs_point_str = obs_point
    else:
        obs_point_str = '[%f,%f,%f]'%(obs_point[0],obs_point[1],obs_point[2])

    columns = ['B_coulomb_x','B_coulomb_y','B_coulomb_z']
    df_name = conf[run+'_derived']+'timeseries/df_coulomb' \
        + '_obs_point=%s_rcut=%f.pkl'%(obs_point_str,rcut)

    dtimes = []
    slice_arrays = []
    for time in list(times):
        outname = conf[run+'_derived'] + 'timeseries/slices/' \
            + 'B_coulomb_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
            + '_obs_point=%s_rcut=%f.npy'%(obs_point_str,rcut)

        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        if os.path.exists(outname):
            slice_arrays.append(np.load(outname))
        else:
            raise ValueError ('no data for time '\
                +str(time)+'. \nrun slice_B_coulomb(...)')

    df = pd.DataFrame(data=slice_arrays, columns = columns,
                        index=dtimes)
    df.to_pickle(df_name)
