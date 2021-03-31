import os
import sys
import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
from units_and_constants import phys
from numba import njit
import datetime
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../swmf/')
import read_swmf_files as rswmf
import cxtransform as cx
from named_var_indexes import index2str, nVarTot, nVarNeeded, \
                            _x                   ,\
                            _y                   ,\
                            _z                   ,\
                            _rho                 ,\
                            _ux                  ,\
                            _uy                  ,\
                            _uz                  ,\
                            _e                   ,\
                            _bx                  ,\
                            _by                  ,\
                            _bz                  ,\
                            _b1x                 ,\
                            _b1y                 ,\
                            _b1z                 ,\
                            _p                   ,\
                            _jx                  ,\
                            _jy                  ,\
                            _jz                  ,\
                            _norm_b              ,\
                            _norm_b1             ,\
                            _norm_j              ,\
                            _div_b               ,\
                            _div_b1              ,\
                            _div_j               ,\
                            _curl_b_x            ,\
                            _curl_b_y            ,\
                            _curl_b_z            ,\
                            _norm_curl_b         ,\
                            _curl_b1_x           ,\
                            _curl_b1_y           ,\
                            _curl_b1_z           ,\
                            _norm_curl_b1        ,\
                            _curl_j_x            ,\
                            _curl_j_y            ,\
                            _curl_j_z            ,\
                            _norm_curl_j         ,\
                            _relative_div_b      ,\
                            _relative_div_b1     ,\
                            _relative_div_j      ,\
                            _FNdel_b             ,\
                            _FNdel_b1            ,\
                            _FNdel_j             ,\
                            _normalized_div_b    ,\
                            _normalized_div_b1   ,\
                            _normalized_div_j    ,\
                            _div_over_curl_b     ,\
                            _div_over_curl_b1    ,\
                            _div_over_curl_j     ,\
                            _jRx                 ,\
                            _jRy                 ,\
                            _jRz                 ,\
                            _norm_jR             ,\
                            _jR_error            ,\
                            _jR_fractional_error ,\
                            _gridspacing


def earth_surface(mlat,mlon,time):
    return cx.MAGtoGSM(np.array([1.,mlat,mlon]),time,'sph','car')


@njit
def B_biotsavart(point, NeededData, rcut=0.):
    print('start')
    unique_epsilons = np.array([  0.0625,
                                  0.1250,
                                  0.2500,
                                  0.5000,
                                  1.0000,
                                  2.0000,
                                  4.0000,
                                  8.0000 ])
    n_eps = unique_epsilons.size

    nBlock, nI, nJ, nK = NeededData.shape[1:]

    ret = np.zeros((n_eps,3), dtype=np.float32)
    for iBlockP in range(nBlock):
        epsilon = NeededData[_x,iBlockP,1,0,0] - NeededData[_x,iBlockP,0,0,0]
        assert(np.argwhere(unique_epsilons==epsilon).shape == (1,1))
        i_eps = np.argwhere(unique_epsilons==epsilon)[0][0]
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):
                    if i == 0 or j == 0 or k == 0 or i == nI-1 or j == nJ-1 or k == nK-1:
                        continue
                    distanceSquared = ( NeededData[_x, iBlockP, i,j,k]**2 \
                                      + NeededData[_y, iBlockP, i,j,k]**2 \
                                      + NeededData[_z, iBlockP, i,j,k]**2 )
                    if distanceSquared < rcut**2:
                        continue

                    partials = np.empty((3,3),dtype=np.float32); partials[:,:]=np.nan
                    integrand = np.empty((3,),dtype=np.float32); integrand[:]=np.nan

                    partials[0,0] = (NeededData[_b1x, iBlockP, i+1, j  , k  ] - NeededData[_b1x, iBlockP, i-1, j  , k  ])/(2*epsilon)
                    partials[0,1] = (NeededData[_b1y, iBlockP, i+1, j  , k  ] - NeededData[_b1y, iBlockP, i-1, j  , k  ])/(2*epsilon)
                    partials[0,2] = (NeededData[_b1z, iBlockP, i+1, j  , k  ] - NeededData[_b1z, iBlockP, i-1, j  , k  ])/(2*epsilon)
                    partials[1,0] = (NeededData[_b1x, iBlockP, i  , j+1, k  ] - NeededData[_b1x, iBlockP, i  , j-1, k  ])/(2*epsilon)
                    partials[1,1] = (NeededData[_b1y, iBlockP, i  , j+1, k  ] - NeededData[_b1y, iBlockP, i  , j-1, k  ])/(2*epsilon)
                    partials[1,2] = (NeededData[_b1z, iBlockP, i  , j+1, k  ] - NeededData[_b1z, iBlockP, i  , j-1, k  ])/(2*epsilon)
                    partials[2,0] = (NeededData[_b1x, iBlockP, i  , j  , k+1] - NeededData[_b1x, iBlockP, i  , j  , k-1])/(2*epsilon)
                    partials[2,1] = (NeededData[_b1y, iBlockP, i  , j  , k+1] - NeededData[_b1y, iBlockP, i  , j  , k-1])/(2*epsilon)
                    partials[2,2] = (NeededData[_b1z, iBlockP, i  , j  , k+1] - NeededData[_b1z, iBlockP, i  , j  , k-1])/(2*epsilon)

                    curl_B1_tens = partials - partials.transpose()
                    curl_B1_x = curl_B1_tens[1,2]
                    curl_B1_y = curl_B1_tens[2,0]
                    curl_B1_z = curl_B1_tens[0,1]

                    r_x = point[0] - NeededData[_x, iBlockP, i, j, k]
                    r_y = point[1] - NeededData[_y, iBlockP, i, j, k]
                    r_z = point[2] - NeededData[_z, iBlockP, i, j, k]
                    r = np.sqrt(r_x**2 + r_y**2 + r_z**2)

                    integrand[2] = curl_B1_x*r_y - curl_B1_y*r_x
                    integrand[0] = curl_B1_y*r_z - curl_B1_z*r_y
                    integrand[1] = curl_B1_z*r_x - curl_B1_x*r_z
                    integrand[:] =  integrand[:]/r**3

                    ret[i_eps,:] = ret[i_eps,:] + integrand

    for i_eps in range(n_eps):
        ret[i_eps,:] = ( (unique_epsilons[i_eps]**3)/(4*np.pi) ) * ret[i_eps,:]

    if False:
        return ret
    else:
        return(np.sum(ret, axis=0))


@njit
def B_coulomb(point, NeededData, rcut=0.):
    print('start')
    unique_epsilons = np.array([  0.0625,
                                  0.1250,
                                  0.2500,
                                  0.5000,
                                  1.0000,
                                  2.0000,
                                  4.0000,
                                  8.0000 ])
    n_eps = unique_epsilons.size

    nBlock, nI, nJ, nK = NeededData.shape[1:]

    ret = np.zeros((n_eps,3), dtype=np.float32)
    for iBlockP in range(nBlock):
        epsilon = NeededData[_x,iBlockP,1,0,0] - NeededData[_x,iBlockP,0,0,0]
        assert(np.argwhere(unique_epsilons==epsilon).shape == (1,1))
        i_eps = np.argwhere(unique_epsilons==epsilon)[0][0]
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):
                    if i == 0 or j == 0 or k == 0 or i == nI-1 or j == nJ-1 or k == nK-1:
                        continue
                    distanceSquared = ( NeededData[_x, iBlockP, i,j,k]**2 \
                                      + NeededData[_y, iBlockP, i,j,k]**2 \
                                      + NeededData[_z, iBlockP, i,j,k]**2 )
                    if distanceSquared < rcut**2:
                        continue

                    partials = np.empty((3,3),dtype=np.float32); partials[:,:]=np.nan
                    integrand = np.empty((3,),dtype=np.float32); integrand[:]=np.nan

                    partials[0,0] = (NeededData[_b1x, iBlockP, i+1, j  , k  ] - NeededData[_b1x, iBlockP, i-1, j  , k  ])/(2*epsilon)
                    partials[0,1] = (NeededData[_b1y, iBlockP, i+1, j  , k  ] - NeededData[_b1y, iBlockP, i-1, j  , k  ])/(2*epsilon)
                    partials[0,2] = (NeededData[_b1z, iBlockP, i+1, j  , k  ] - NeededData[_b1z, iBlockP, i-1, j  , k  ])/(2*epsilon)
                    partials[1,0] = (NeededData[_b1x, iBlockP, i  , j+1, k  ] - NeededData[_b1x, iBlockP, i  , j-1, k  ])/(2*epsilon)
                    partials[1,1] = (NeededData[_b1y, iBlockP, i  , j+1, k  ] - NeededData[_b1y, iBlockP, i  , j-1, k  ])/(2*epsilon)
                    partials[1,2] = (NeededData[_b1z, iBlockP, i  , j+1, k  ] - NeededData[_b1z, iBlockP, i  , j-1, k  ])/(2*epsilon)
                    partials[2,0] = (NeededData[_b1x, iBlockP, i  , j  , k+1] - NeededData[_b1x, iBlockP, i  , j  , k-1])/(2*epsilon)
                    partials[2,1] = (NeededData[_b1y, iBlockP, i  , j  , k+1] - NeededData[_b1y, iBlockP, i  , j  , k-1])/(2*epsilon)
                    partials[2,2] = (NeededData[_b1z, iBlockP, i  , j  , k+1] - NeededData[_b1z, iBlockP, i  , j  , k-1])/(2*epsilon)

                    div_B1 = partials[0,0] + partials[1,1] + partials[2,2]

                    r_x = point[0] - NeededData[_x, iBlockP, i, j, k]
                    r_y = point[1] - NeededData[_y, iBlockP, i, j, k]
                    r_z = point[2] - NeededData[_z, iBlockP, i, j, k]
                    r = np.sqrt(r_x**2 + r_y**2 + r_z**2)

                    integrand[0] = div_B1*r_x
                    integrand[1] = div_B1*r_y
                    integrand[2] = div_B1*r_z
                    integrand[:] =  integrand[:]/r**3

                    ret[i_eps,:] = ret[i_eps,:] + integrand

    for i_eps in range(n_eps):
        ret[i_eps,:] = ( (unique_epsilons[i_eps]**3)/(4*np.pi) ) * ret[i_eps,:]

    if False:
        return ret
    else:
        return(np.sum(ret, axis=0))


def dst(run, para=True, n_times=50, debug=False, rcut=0.):
    times = list(util.get_available_slices(run)[1])
    if n_times is not None:
        times = times[:10*n_times:10]

    def dst_wrap(time):
        NeededData = rswmf.get_needed_array(run,time)
        Bbs = B_biotsavart( np.zeros(3), NeededData, rcut=rcut)
        Bcl = B_coulomb(    np.zeros(3), NeededData, rcut=rcut)
        return (np.linalg.norm(Bbs),np.linalg.norm(Bcl))

    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is not None and num_cores > len(times):
            num_cores = len(times)
        print('Parallel processing {0:d} time slices using {1:d} cores'\
              .format(len(times), num_cores))
        dsts = Parallel(n_jobs=num_cores)(\
                        delayed(dst_wrap)(time) for time in times)
    else:
        dsts = []
        for time in times:
            dsts.append(dst_wrap(time))

    dtimes = []
    for time in times:
        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))

    print(dsts)
    df_dst = pd.DataFrame(data=dsts, index=dtimes, columns=['dst_biotsavart', 'dst_coulomb'])
    df_dst.to_pickle('df_dst.pkl') 

if __name__=='__main__':
    run = 'DIPTSUR2'
    time = (2019,9,2,4,10,0,0)

    #point = earth_surface(0,0,time)
    #point = np.zeros(3)
    #
    #Bbs = B_biotsavart(point,
    #                    rswmf.get_needed_array(run,time),
    #                    rcut=util.get_rCurrents(run))
    #print(Bbs)
    #
    #Bcl = B_coulomb(point,
    #                    rswmf.get_needed_array(run,time),
    #                    rcut=util.get_rCurrents(run))
    #print(Bcl)

    dst(run, para=True, rcut=util.get_rCurrents(run))
