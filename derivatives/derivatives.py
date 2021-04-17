import numpy as np
from numba import njit
from named_var_indexes import _x,_b1x,_b1y,_b1z,_bx,_by,_bz,_jx,_jy,_jz,_ux,_uy,_uz

@njit
def get_partials(DataArray, vvar, iBlockP,i,j,k):
    partials = np.empty((3,3),dtype=np.float32); partials[:,:]=np.nan

    nBlock, nI, nJ, nK = DataArray.shape[1:]
    epsilon = DataArray[_x,iBlockP,1,0,0] - DataArray[_x,iBlockP,0,0,0]

    if vvar=='u':
        _fx = _ux
        _fy = _uy
        _fz = _uz
    if vvar=='b1':
        _fx = _b1x
        _fy = _b1y
        _fz = _b1z
    if vvar=='b':
        _fx = _bx
        _fy = _by
        _fz = _bz
    if vvar=='j':
        _fx = _jx
        _fy = _jy
        _fz = _jz

    if i == 0 or j == 0 or k == 0 or i == nI-1 or j == nJ-1 or k == nK-1:
        return partials

    partials[0,0] = (DataArray[_fx, iBlockP, i+1, j  , k  ] - DataArray[_fx, iBlockP, i-1, j  , k  ])/(2*epsilon)
    partials[0,1] = (DataArray[_fy, iBlockP, i+1, j  , k  ] - DataArray[_fy, iBlockP, i-1, j  , k  ])/(2*epsilon)
    partials[0,2] = (DataArray[_fz, iBlockP, i+1, j  , k  ] - DataArray[_fz, iBlockP, i-1, j  , k  ])/(2*epsilon)
    partials[1,0] = (DataArray[_fx, iBlockP, i  , j+1, k  ] - DataArray[_fx, iBlockP, i  , j-1, k  ])/(2*epsilon)
    partials[1,1] = (DataArray[_fy, iBlockP, i  , j+1, k  ] - DataArray[_fy, iBlockP, i  , j-1, k  ])/(2*epsilon)
    partials[1,2] = (DataArray[_fz, iBlockP, i  , j+1, k  ] - DataArray[_fz, iBlockP, i  , j-1, k  ])/(2*epsilon)
    partials[2,0] = (DataArray[_fx, iBlockP, i  , j  , k+1] - DataArray[_fx, iBlockP, i  , j  , k-1])/(2*epsilon)
    partials[2,1] = (DataArray[_fy, iBlockP, i  , j  , k+1] - DataArray[_fy, iBlockP, i  , j  , k-1])/(2*epsilon)
    partials[2,2] = (DataArray[_fz, iBlockP, i  , j  , k+1] - DataArray[_fz, iBlockP, i  , j  , k-1])/(2*epsilon)
    return partials
