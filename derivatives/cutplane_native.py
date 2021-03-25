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


@njit(error_model='numpy')
def xzplane(NeededData):
    #16384==128**2
    print('starting')
    planedata = np.empty((4,16384), dtype=np.float32); planedata[:,:]=np.nan
    nBlock, nI, nJ, nK = NeededData.shape[1:]
    print('check1')
    counter = 0
    print('check2')

    for iBlockP in range(nBlock):
        epsilon = NeededData[_x,iBlockP,1,0,0] - NeededData[_x,iBlockP,0,0,0]
        if epsilon != 0.0625:
            continue
        for j in range(nJ):
            if 0.09375 != NeededData[_y, iBlockP, 0,j,0]:
                continue
            for i in range(nI):
                for k in range(nK):
                    planedata[0,counter] = NeededData[_x,iBlockP,i,j,k]
                    planedata[1,counter] = NeededData[_z,iBlockP,i,j,k]
                    if i != 0 and j != 0 and k != 0 and i != nI-1 and j != nJ-1 and k != nK-1:
                        planedata[2,counter] = (NeededData[_b1x, iBlockP, i+1, j  , k  ] - NeededData[_b1x, iBlockP, i-1, j  , k  ])/(2*epsilon) \
                                             + (NeededData[_b1y, iBlockP, i  , j+1, k  ] - NeededData[_b1y, iBlockP, i  , j-1, k  ])/(2*epsilon) \
                                             + (NeededData[_b1z, iBlockP, i  , j  , k+1] - NeededData[_b1z, iBlockP, i  , j  , k-1])/(2*epsilon)

                    planedata[3,counter] = np.sqrt( NeededData[_b1x,iBlockP,i,j,k]**2 \
                                                  + NeededData[_b1y,iBlockP,i,j,k]**2 \
                                                  + NeededData[_b1z,iBlockP,i,j,k]**2 )

                    counter += 1

    print(counter)
    return planedata

def xzplane_wrap(run, time):
    data, dTREE, ind, block2node, node2block = rswmf.read_all(util.time2CDFfilename(run,time)[:-8])
    nBlock, nI, nJ, nK = block2node.size, dTREE['nI'], dTREE['nJ'], dTREE['nK']
    NeededData = np.empty((nVarNeeded, nBlock, nI, nJ, nK), dtype=np.float32); NeededData[:,:,:,:,:] = np.nan
    for index in range(nVarNeeded):
        NeededData[index,:,:,:,:] = data[index2str[index]][ind].reshape((nBlock, nI, nJ, nK))
    del data, dTREE, ind, block2node, node2block

    direct = '/home/gary/magnetosphere/images/'+run+'/cutplane_native/'
    if not os.path.exists(direct): os.makedirs(direct)
    fname = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_x_z_divB1_normB1.npy'%util.tpad(time, length=6)
    planedata = xzplane(NeededData)
    print('saving '+fname)
    np.save(fname,planedata)


if __name__ == '__main__':
    ##################
    run = 'DIPTSUR2'
    debug = True

    if run == 'DIPTSUR2':
        time = (2019,9,2,6,30,0,0)
        #time = (2019,9,2,4,10,0,0)
    if run == 'IMP10_RUN_SAMPLE':
        time = (2019,9,2,7,0,0,0)
    if run == 'TESTANALYTIC':
        time = (2000,1,1,0,10,0,0)
    if run == 'LUHMANN1979':
        time = (2000,1,1,0,0,0,0)
    ##################

    xzplane_wrap(run, (2019,9,2,6,30,0,0))
    xzplane_wrap(run, (2019,9,2,4,10,0,0))

