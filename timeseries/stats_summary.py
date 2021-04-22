#python -c "import stats_summary as s; s.slice_stats_summary('DIPTSUR2',(2019,9,2,4,10,0,0))"
#python -c "import stats_summary as s; s.stitch_stats_summary('DIPTSUR2',[(2019,9,2,4,10,0,0)])"
import numpy as np
from numba import njit
import pandas as pd
import datetime
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
from swmf_file_reader.read_swmf_files import read_all
from derivatives import get_partials
from units_and_constants import phys
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
_count = 0
_mean  = 1
_sndMo = 2
_std   = 3
_min   = 4
_max   = 5

@njit(error_model='numpy')
def _jit_stats_summary(DataArray, rcut):
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

    summary_arr = np.empty((nVarTot,n_eps,6), dtype=np.float32)

    summary_arr[:,:,_count] = 0.
    summary_arr[:,:,_mean ] = 0.
    summary_arr[:,:,_sndMo] = 0.
    summary_arr[:,:,_std  ] = 0.
    summary_arr[:,:,_min  ] = np.inf
    summary_arr[:,:,_max  ] = -np.inf

    store = np.empty((nVarTot,), dtype=np.float32)

    for iBlockP in range(nBlock):
        epsilon = DataArray[_x,iBlockP,1,0,0] - DataArray[_x,iBlockP,0,0,0]
        assert(np.argwhere(unique_epsilons==epsilon).shape == (1,1))
        i_eps = np.argwhere(unique_epsilons==epsilon)[0][0]
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):
                    store[:]=np.nan
                    distanceSquared = ( DataArray[_x, iBlockP, i,j,k]**2 \
                                      + DataArray[_y, iBlockP, i,j,k]**2 \
                                      + DataArray[_z, iBlockP, i,j,k]**2 )
                    if distanceSquared <= rcut**2:
                        continue

                    store[_gridspacing] = epsilon

                    for i_var in range(nVarNeeded):
                        store[i_var] = DataArray[i_var, iBlockP, i,j,k]

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

                        store[_norm_f] = np.sqrt( store[_fx]**2 \
                                                + store[_fy]**2 \
                                                + store[_fz]**2  )

                        partials = get_partials(DataArray, vvar, iBlockP,i,j,k)

                        store[_div_f] = partials[0,0]+partials[1,1]+partials[2,2]

                        curl_tens = partials - partials.transpose()
                        store[_curl_f_x] = curl_tens[1,2]
                        store[_curl_f_y] = curl_tens[2,0]
                        store[_curl_f_z] = curl_tens[0,1]
                        store[_norm_curl_f] = np.sqrt( store[_curl_f_x]**2 \
                                                     + store[_curl_f_y]**2 \
                                                     + store[_curl_f_z]**2  )

                        store[_FNdel_f] = np.sqrt(np.sum(partials**2))

                        store[_relative_div_f] = store[_div_f]/store[_norm_f]
                        store[_normalized_div_f] = store[_div_f]/store[_FNdel_f]
                        store[_div_over_curl_f] = store[_div_f]/store[_norm_curl_f]

                    store[_jRx] = (1./unitmu0)*store[_curl_b1_x]
                    store[_jRy] = (1./unitmu0)*store[_curl_b1_y]
                    store[_jRz] = (1./unitmu0)*store[_curl_b1_z]
                    store[_norm_jR] = np.sqrt( store[_jRx]**2 \
                                             + store[_jRy]**2 \
                                             + store[_jRz]**2 )
                    store[_jR_error] = np.sqrt( (store[_jx]-store[_jRx])**2 \
                                              + (store[_jy]-store[_jRy])**2 \
                                              + (store[_jz]-store[_jRz])**2 )
                    store[_jR_fractional_error] = store[_jR_error]/store[_norm_jR]

                    for i_var in range(nVarTot):
                        if not np.isnan(store[i_var]):
                            summary_arr[i_var,i_eps,_count] = summary_arr[i_var,i_eps,_count] + 1
                            summary_arr[i_var,i_eps,_mean ] = summary_arr[i_var,i_eps,_mean ] + store[i_var]
                            summary_arr[i_var,i_eps,_sndMo] = summary_arr[i_var,i_eps,_sndMo] + store[i_var]**2
                            if store[i_var] < summary_arr[i_var,i_eps,_min]:
                                summary_arr[i_var,i_eps,_min] = store[i_var]
                            if store[i_var] > summary_arr[i_var,i_eps,_max]:
                                summary_arr[i_var,i_eps,_max] = store[i_var]

    summary_arr[:,:,_mean ] = summary_arr[:,:,_mean ]/summary_arr[:,:,_count]
    summary_arr[:,:,_sndMo] = summary_arr[:,:,_sndMo]/summary_arr[:,:,_count]
    summary_arr[:,:,_std  ] = np.sqrt(summary_arr[:,:,_sndMo] - summary_arr[:,:,_mean]**2)#!! could suffer catastrophic cancelation
    return summary_arr


def slice_stats_summary(run, time, rcut=None, cache=None):
    if cache is None:
        cache = read_all(util.time2CDFfilename(run,time)[:-8])

    if rcut is None:
        rcut = util.get_rCurrents(run)

    outname = conf[run+'_derived'] + 'timeseries/slices/' \
        + 'stats_summary_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
        + '_rcut=%f.npy'%(rcut)
    summary_arr = _jit_stats_summary(cache['DataArray'], rcut)
    np.save(outname, summary_arr)


def stitch_stats_summary(run, times, rcut=None):
    if rcut is None:
        rcut = util.get_rCurrents(run)

    epsilons = [ 0.0625,
                 0.1250,
                 0.2500,
                 0.5000,
                 1.0000,
                 2.0000,
                 4.0000,
                 8.0000 ]

    direct = conf[run+'_derived'] + 'timeseries/'

    list_summaries = []
    dtimes = []
    for time in list(times):
        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        outname = conf[run+'_derived'] + 'timeseries/slices/' \
            + 'stats_summary_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
            + '_rcut=%f.npy'%(rcut)
        list_summaries.append(np.load(outname))

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

        summaryDF.to_pickle(direct+'stats_summary_%s.pkl'%(index2str[i_var]))
