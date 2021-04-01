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


unitmu0 = np.float32( phys['mu0']*(phys['muA']/phys['m']**2) )

@njit(error_model='numpy')
def get_summary(NeededData, rcut):
    unique_epsilons = np.array([  0.0625,
                                  0.1250,
                                  0.2500,
                                  0.5000,
                                  1.0000,
                                  2.0000,
                                  4.0000,
                                  8.0000 ])
    n_eps = unique_epsilons.size

    _count = 0
    _mean  = 1
    _sndMo = 2
    _std   = 3
    _min   = 4
    _max   = 5

    nBlock, nI, nJ, nK = NeededData.shape[1:]

    summary_arr = np.empty((nVarTot,n_eps,6), dtype=np.float32)

    summary_arr[:,:,_count] = 0.
    summary_arr[:,:,_mean ] = 0.
    summary_arr[:,:,_sndMo] = 0.
    summary_arr[:,:,_std  ] = 0.
    summary_arr[:,:,_min  ] = np.inf
    summary_arr[:,:,_max  ] = -np.inf

    cache = np.empty((nVarTot,), dtype=np.float32)
    partials = np.empty((3,3), dtype=np.float32)

    for iBlockP in range(nBlock):
        epsilon = NeededData[_x,iBlockP,1,0,0] - NeededData[_x,iBlockP,0,0,0]
        assert(np.argwhere(unique_epsilons==epsilon).shape == (1,1))
        i_eps = np.argwhere(unique_epsilons==epsilon)[0][0]
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):
                    cache[:]=np.nan
                    distanceSquared = ( NeededData[_x, iBlockP, i,j,k]**2 \
                                      + NeededData[_y, iBlockP, i,j,k]**2 \
                                      + NeededData[_z, iBlockP, i,j,k]**2 )
                    if distanceSquared <= rcut**2:
                        continue

                    cache[_gridspacing] = epsilon

                    for i_var in range(nVarNeeded):
                        cache[i_var]   = NeededData[i_var, iBlockP, i,j,k]

                    for vvar in ['b','b1','j']:
                        if vvar == 'b':
                            _fx               = _bx              
                            _fy               = _by              
                            _fz               = _bz              
                            _norm_f           = _norm_b          
                            _div_f            = _div_b           
                            _curl_f_x         = _curl_b_x        
                            _curl_f_y         = _curl_b_y        
                            _curl_f_z         = _curl_b_z        
                            _norm_curl_f      = _norm_curl_b     
                            _FNdel_f          = _FNdel_b         
                            _relative_div_f   = _relative_div_b  
                            _normalized_div_f = _normalized_div_b
                            _div_over_curl_f  = _div_over_curl_b 
                        if vvar == 'b1':
                            _fx               = _b1x              
                            _fy               = _b1y              
                            _fz               = _b1z              
                            _norm_f           = _norm_b1          
                            _div_f            = _div_b1           
                            _curl_f_x         = _curl_b1_x        
                            _curl_f_y         = _curl_b1_y        
                            _curl_f_z         = _curl_b1_z        
                            _norm_curl_f      = _norm_curl_b1     
                            _FNdel_f          = _FNdel_b1         
                            _relative_div_f   = _relative_div_b1  
                            _normalized_div_f = _normalized_div_b1
                            _div_over_curl_f  = _div_over_curl_b1 
                        if vvar == 'j':
                            _fx               = _jx              
                            _fy               = _jy              
                            _fz               = _jz              
                            _norm_f           = _norm_j          
                            _div_f            = _div_j           
                            _curl_f_x         = _curl_j_x        
                            _curl_f_y         = _curl_j_y        
                            _curl_f_z         = _curl_j_z        
                            _norm_curl_f      = _norm_curl_j     
                            _FNdel_f          = _FNdel_j         
                            _relative_div_f   = _relative_div_j  
                            _normalized_div_f = _normalized_div_j
                            _div_over_curl_f  = _div_over_curl_j 

                        cache[_norm_f] = np.sqrt( cache[_fx]**2 \
                                                + cache[_fy]**2 \
                                                + cache[_fz]**2  )

                        partials[:,:] = np.nan
                        if i != 0 and j != 0 and k != 0 and i != nI-1 and j != nJ-1 and k != nK-1:
                            partials[0,0] = (NeededData[_fx, iBlockP, i+1, j  , k  ] - NeededData[_fx, iBlockP, i-1, j  , k  ])/(2*epsilon)
                            partials[0,1] = (NeededData[_fy, iBlockP, i+1, j  , k  ] - NeededData[_fy, iBlockP, i-1, j  , k  ])/(2*epsilon)
                            partials[0,2] = (NeededData[_fz, iBlockP, i+1, j  , k  ] - NeededData[_fz, iBlockP, i-1, j  , k  ])/(2*epsilon)
                            partials[1,0] = (NeededData[_fx, iBlockP, i  , j+1, k  ] - NeededData[_fx, iBlockP, i  , j-1, k  ])/(2*epsilon)
                            partials[1,1] = (NeededData[_fy, iBlockP, i  , j+1, k  ] - NeededData[_fy, iBlockP, i  , j-1, k  ])/(2*epsilon)
                            partials[1,2] = (NeededData[_fz, iBlockP, i  , j+1, k  ] - NeededData[_fz, iBlockP, i  , j-1, k  ])/(2*epsilon)
                            partials[2,0] = (NeededData[_fx, iBlockP, i  , j  , k+1] - NeededData[_fx, iBlockP, i  , j  , k-1])/(2*epsilon)
                            partials[2,1] = (NeededData[_fy, iBlockP, i  , j  , k+1] - NeededData[_fy, iBlockP, i  , j  , k-1])/(2*epsilon)
                            partials[2,2] = (NeededData[_fz, iBlockP, i  , j  , k+1] - NeededData[_fz, iBlockP, i  , j  , k-1])/(2*epsilon)

                        cache[_div_f] = partials[0,0]+partials[1,1]+partials[2,2]

                        curl_tens = partials - partials.transpose()
                        cache[_curl_f_x] = curl_tens[1,2]
                        cache[_curl_f_y] = curl_tens[2,0]
                        cache[_curl_f_z] = curl_tens[0,1]
                        cache[_norm_curl_f] = np.sqrt( cache[_curl_f_x]**2 \
                                                     + cache[_curl_f_y]**2 \
                                                     + cache[_curl_f_z]**2  )

                        cache[_FNdel_f] = np.sqrt(np.sum(partials**2))

                        cache[_relative_div_f] = cache[_div_f]/cache[_norm_f]
                        cache[_normalized_div_f] = cache[_div_f]/cache[_FNdel_f]
                        cache[_div_over_curl_f] = cache[_div_f]/cache[_norm_curl_f]

                    cache[_jRx] = (1./unitmu0)*cache[_curl_b1_x]
                    cache[_jRy] = (1./unitmu0)*cache[_curl_b1_y]
                    cache[_jRz] = (1./unitmu0)*cache[_curl_b1_z]
                    cache[_norm_jR] = np.sqrt( cache[_jRx]**2 \
                                             + cache[_jRy]**2 \
                                             + cache[_jRz]**2 )
                    cache[_jR_error] = np.sqrt( (cache[_jx]-cache[_jRx])**2 \
                                              + (cache[_jy]-cache[_jRy])**2 \
                                              + (cache[_jz]-cache[_jRz])**2 )
                    cache[_jR_fractional_error] = cache[_jR_error]/cache[_norm_jR]

                    for i_var in range(nVarTot):
                        if not np.isnan(cache[i_var]):
                            summary_arr[i_var,i_eps,_count] = summary_arr[i_var,i_eps,_count] + 1
                            summary_arr[i_var,i_eps,_mean ] = summary_arr[i_var,i_eps,_mean ] + cache[i_var]
                            summary_arr[i_var,i_eps,_sndMo] = summary_arr[i_var,i_eps,_sndMo] + cache[i_var]**2
                            if cache[i_var] < summary_arr[i_var,i_eps,_min]:
                                summary_arr[i_var,i_eps,_min] = cache[i_var]
                            if cache[i_var] > summary_arr[i_var,i_eps,_max]:
                                summary_arr[i_var,i_eps,_max] = cache[i_var]

    summary_arr[:,:,_mean ] = summary_arr[:,:,_mean ]/summary_arr[:,:,_count]
    summary_arr[:,:,_sndMo] = summary_arr[:,:,_sndMo]/summary_arr[:,:,_count]
    summary_arr[:,:,_std  ] = np.sqrt(summary_arr[:,:,_sndMo] - summary_arr[:,:,_mean]**2)

    return summary_arr


def summarize(run, time, debug=False, rcut=2.):
    summary_arr = get_summary(rswmf.get_needed_array(run,time), rcut)

    direct = conf[run+'_derived'] + 'derivatives/native_grid/'
    if not os.path.exists(direct): os.makedirs(direct)
    fname_summary = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_summary.npy'%util.tpad(time, length=6)

    print('saving '+fname_summary)
    np.save(fname_summary,summary_arr)

def summarize_all(run, para=True, n_times=None, debug=False, rcut=2.):
    times = list(util.get_available_slices(run)[1])
    if n_times is not None:
        times = times[:n_times]

    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is not None and num_cores > len(times):
            num_cores = len(times)
        print('Parallel processing {0:d} time slices using {1:d} cores'\
              .format(len(times), num_cores))
        Parallel(n_jobs=num_cores)(\
                delayed(summarize)(run, time, debug=debug, rcut=rcut) for time in times)
    else:
        for time in times:
            summarize(run, time, debug=debug, rcut=rcut)

def stitch_together(run, n_times=None):
    epsilons = [ 0.0625,
                 0.1250,
                 0.2500,
                 0.5000,
                 1.0000,
                 2.0000,
                 4.0000,
                 8.0000 ]
    _count = 0
    _mean  = 1
    _sndMo = 2
    _std   = 3
    _min   = 4
    _max   = 5 #!!! warning: duplicated

    times = list(util.get_available_slices(run)[1])
    if n_times is not None:
        times = times[:n_times]

    direct = conf[run+'_derived'] + 'derivatives/native_grid/'
    list_summaries = []
    dtimes = []
    for time in times:
        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        fname_summary = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_summary.npy'%util.tpad(time, length=6)
        list_summaries.append(np.load(fname_summary))

    multiind = pd.MultiIndex.from_product([epsilons, ['count', 'mean', 'sndMo', 'std', 'min', 'max']])
    for i_var in range(nVarTot):
        summaryDF = pd.DataFrame(index=dtimes, columns=multiind)

        for i_eps in range(len(epsilons)):
            for i_t in range(len(times)):
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'count')] = list_summaries[i_t][i_var,i_eps, _count]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'mean' )] = list_summaries[i_t][i_var,i_eps, _mean ]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'sndMo')] = list_summaries[i_t][i_var,i_eps, _sndMo]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'std'  )] = list_summaries[i_t][i_var,i_eps, _std  ]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'min'  )] = list_summaries[i_t][i_var,i_eps, _min  ]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'max'  )] = list_summaries[i_t][i_var,i_eps, _max  ]

        summaryDF.to_pickle(direct+'summaryDF_%s.pkl'%(index2str[i_var]))

if __name__ == '__main__':
    ##################
    run = 'DIPTSUR2'
    debug = True

    if run == 'DIPTSUR2':
        #time = (2019,9,2,6,30,0,0)
        time = (2019,9,2,4,10,0,0)
    if run == 'IMP10_RUN_SAMPLE':
        time = (2019,9,2,7,0,0,0)
    if run == 'TESTANALYTIC':
        time = (2000,1,1,0,10,0,0)
    if run == 'LUHMANN1979':
        time = (2000,1,1,0,0,0,0)
    ##################

    #main(run, time, debug=debug)

    n_times = None
    para = True
    summarize_all(run,para=para,n_times=n_times,debug=debug)
    stitch_together(run, n_times=n_times)
