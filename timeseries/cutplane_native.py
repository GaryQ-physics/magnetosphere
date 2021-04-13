import os
import sys
import numpy as np
import pandas as pd
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
from units_and_constants import phys
from numba import njit
import datetime
from read_swmf_files2 import read_all
from derivatives import get_partials

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
def _jit_xzplane(DataArray, resolution, expectedNpts):
    planedata = np.empty((5,expectedNpts), dtype=np.float32); planedata[:,:]=np.nan
    nBlock, nI, nJ, nK = DataArray.shape[1:]
    counter = 0

    for iBlockP in range(nBlock):
        epsilon = DataArray[_x,iBlockP,1,0,0] - DataArray[_x,iBlockP,0,0,0]
        if epsilon != resolution:
            continue
        for j in range(nJ):
            if resolution*3./2. != DataArray[_y, iBlockP, 0,j,0]:
                continue
            for i in range(nI):
                for k in range(nK):
                    if counter >= expectedNpts:
                        print('WARNING: exceeded expectedNpts')
                        break
                    planedata[0,counter] = DataArray[_x,iBlockP,i,j,k]
                    planedata[1,counter] = DataArray[_z,iBlockP,i,j,k]

                    partials = get_partials(DataArray, 'b1', iBlockP,i,j,k)
                   # if i != 0 and j != 0 and k != 0 and i != nI-1 and j != nJ-1 and k != nK-1:
                   #     planedata[2,counter] = (NeededData[_b1x, iBlockP, i+1, j  , k  ] - NeededData[_b1x, iBlockP, i-1, j  , k  ])/(2*epsilon) \
                   #                          + (NeededData[_b1y, iBlockP, i  , j+1, k  ] - NeededData[_b1y, iBlockP, i  , j-1, k  ])/(2*epsilon) \
                   #                          + (NeededData[_b1z, iBlockP, i  , j  , k+1] - NeededData[_b1z, iBlockP, i  , j  , k-1])/(2*epsilon)
                    planedata[2,counter] = partials[0,0]+partials[1,1]+partials[2,2]

                    planedata[3,counter] = np.sqrt( DataArray[_b1x,iBlockP,i,j,k]**2 \
                                                  + DataArray[_b1y,iBlockP,i,j,k]**2 \
                                                  + DataArray[_b1z,iBlockP,i,j,k]**2 )

                    curl_B1_tens = partials - partials.transpose()
                    curl_B1_x = curl_B1_tens[1,2]
                    curl_B1_y = curl_B1_tens[2,0]
                    curl_B1_z = curl_B1_tens[0,1]

                    planedata[4,counter] = np.sqrt( curl_B1_x**2 \
                                                  + curl_B1_y**2 \
                                                  + curl_B1_z**2 )

                    counter += 1
    return planedata

def slice_xzplane(run, time, rcut=None, png=True, cache=None):
    if cache is None:
        cache = read_all(util.time2CDFfilename(run,time)[:-8])

    if rcut is None:
        rcut = util.get_rCurrents(run)

    inner_planedata = _jit_xzplane(cache['DataArray'], 0.0625, 15360)
    outer_planedata = _jit_xzplane(cache['DataArray'], 0.125 , 17664)

    if png:
        import matplotlib.pyplot as plt
        #https://stackoverflow.com/questions/17201172/a-logarithmic-colorbar-in-matplotlib-scatter-plot
        import matplotlib.colors as mplc

        ## TODO:
        #png_b1 = conf[run+'_derived'] + 'timeseries/slices/' \
        #    + 'xzplane_normb1_%.2d%.2d%.2dT%.2d%.2d%.2d.png'%util.tpad(time, length=6)

        png_derivs = conf[run+'_derived'] + 'timeseries/slices/' \
            + 'xzplane_derivsb1_%.2d%.2d%.2dT%.2d%.2d%.2d.png'%util.tpad(time, length=6)

        title = '%s at %.2d%.2d%.2dT%.2d%.2d%.2d.png'%(run,*util.tpad(time, length=6))

        x =          np.concatenate([inner_planedata[0,:],outer_planedata[0,:]])# + 1e-25
        z =          np.concatenate([inner_planedata[1,:],outer_planedata[1,:]])# + 1e-25
        divB1 =      np.concatenate([inner_planedata[2,:],outer_planedata[2,:]])# + 1e-25
        normB1 =     np.concatenate([inner_planedata[3,:],outer_planedata[3,:]])# + 1e-25
        normcurlB1 = np.concatenate([inner_planedata[4,:],outer_planedata[4,:]])# + 1e-25

        tr =  x**2 + z**2 <= rcut**2

        x[tr] = np.nan
        z[tr] = np.nan
        divB1[tr] = np.nan
        divB1[tr] = np.nan

        fig, axes = plt.subplots(figsize=(32,6), nrows=1, ncols=2,dpi=300)

        sc0=axes[0].scatter(x,z,c=np.abs(divB1), norm=mplc.LogNorm(vmin=1e-5, vmax=1e+4))
        #axes[0].set_colorbar(label='|div_b1| [$nT R_E^{-1}$]')
        axes[0].set_title('|div_b1|')
        axes[0].set_xlabel('x [R_E]')
        axes[0].set_ylabel('z [R_E]')
        fig.colorbar(sc0, ax=axes[0], label='[$nT R_E^{-1}$]')

        sc1=axes[1].scatter(x,z,c=normcurlB1, norm=mplc.LogNorm(vmin=1e-5, vmax=1e+4))
        #axes[1].set_colorbar(label='norm_curl_b1 [$nT R_E^{-1}$]')
        axes[1].set_title('norm_curl_b1')
        axes[1].set_xlabel('x [R_E]')
        axes[1].set_ylabel('z [R_E]')
        fig.colorbar(sc1, ax=axes[1], label='[$nT R_E^{-1}$]')

        fig.suptitle(title)
        fig.savefig(png_derivs)
        plt.clf()
        del plt, fig, axes

    else:
        inner_outname = conf[run+'_derived'] + 'timeseries/slices/' \
            + 'y=3_32_xzplane_%.2d%.2d%.2dT%.2d%.2d%.2d_x_z_divB1_normB1_normcurlB1.npy'%util.tpad(time, length=6)
        outer_outname = conf[run+'_derived'] + 'timeseries/slices/' \
            + 'y=3_16_xzplane_%.2d%.2d%.2dT%.2d%.2d%.2d_x_z_divB1_normB1_normcurlB1.npy'%util.tpad(time, length=6)
        np.save(inner_outname,inner_planedata)
        np.save(outer_outname,outer_planedata)


def stitch_xzplane(run, times, rcut=None, png=True):
    pass


if __name__ == '__main__':
    ##################
    run = 'DIPTSUR2'
    ##################

    slice_xzplane(run, (2019,9,2,6,30,0,0))
    slice_xzplane(run, (2019,9,2,4,10,0,0))
