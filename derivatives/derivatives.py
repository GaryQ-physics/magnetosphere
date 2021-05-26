import numpy as np
from numba import njit
from swmf_file_reader.named_var_indexes import _x,_b1x,_b1y,_b1z,_bx,_by,_bz,_jx,_jy,_jz,_ux,_uy,_uz

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

    if i == 0:
        partials[0,0] = (DataArray[_fx, iBlockP, 1, j  , k  ] - DataArray[_fx, iBlockP, 0, j  , k  ])/(epsilon)
        partials[0,1] = (DataArray[_fy, iBlockP, 1, j  , k  ] - DataArray[_fy, iBlockP, 0, j  , k  ])/(epsilon)
        partials[0,2] = (DataArray[_fz, iBlockP, 1, j  , k  ] - DataArray[_fz, iBlockP, 0, j  , k  ])/(epsilon)
    elif i == nI-1:
        partials[0,0] = (DataArray[_fx, iBlockP, nI-1, j  , k  ] - DataArray[_fx, iBlockP, nI-2, j  , k  ])/(epsilon)
        partials[0,1] = (DataArray[_fy, iBlockP, nI-1, j  , k  ] - DataArray[_fy, iBlockP, nI-2, j  , k  ])/(epsilon)
        partials[0,2] = (DataArray[_fz, iBlockP, nI-1, j  , k  ] - DataArray[_fz, iBlockP, nI-2, j  , k  ])/(epsilon)
    else:
        partials[0,0] = (DataArray[_fx, iBlockP, i+1, j  , k  ] - DataArray[_fx, iBlockP, i-1, j  , k  ])/(2*epsilon)
        partials[0,1] = (DataArray[_fy, iBlockP, i+1, j  , k  ] - DataArray[_fy, iBlockP, i-1, j  , k  ])/(2*epsilon)
        partials[0,2] = (DataArray[_fz, iBlockP, i+1, j  , k  ] - DataArray[_fz, iBlockP, i-1, j  , k  ])/(2*epsilon)

    if j == 0:
        partials[1,0] = (DataArray[_fx, iBlockP, i  , 1, k  ] - DataArray[_fx, iBlockP, i  , 0, k  ])/(epsilon)
        partials[1,1] = (DataArray[_fy, iBlockP, i  , 1, k  ] - DataArray[_fy, iBlockP, i  , 0, k  ])/(epsilon)
        partials[1,2] = (DataArray[_fz, iBlockP, i  , 1, k  ] - DataArray[_fz, iBlockP, i  , 0, k  ])/(epsilon)
    elif j == nJ-1:
        partials[1,0] = (DataArray[_fx, iBlockP, i  , nJ-1, k  ] - DataArray[_fx, iBlockP, i  , nJ-2, k  ])/(epsilon)
        partials[1,1] = (DataArray[_fy, iBlockP, i  , nJ-1, k  ] - DataArray[_fy, iBlockP, i  , nJ-2, k  ])/(epsilon)
        partials[1,2] = (DataArray[_fz, iBlockP, i  , nJ-1, k  ] - DataArray[_fz, iBlockP, i  , nJ-2, k  ])/(epsilon)
    else:
        partials[1,0] = (DataArray[_fx, iBlockP, i  , j+1, k  ] - DataArray[_fx, iBlockP, i  , j-1, k  ])/(2*epsilon)
        partials[1,1] = (DataArray[_fy, iBlockP, i  , j+1, k  ] - DataArray[_fy, iBlockP, i  , j-1, k  ])/(2*epsilon)
        partials[1,2] = (DataArray[_fz, iBlockP, i  , j+1, k  ] - DataArray[_fz, iBlockP, i  , j-1, k  ])/(2*epsilon)

    if k == 0:
        partials[2,0] = (DataArray[_fx, iBlockP, i  , j  , 1] - DataArray[_fx, iBlockP, i  , j  , 0])/(epsilon)
        partials[2,1] = (DataArray[_fy, iBlockP, i  , j  , 1] - DataArray[_fy, iBlockP, i  , j  , 0])/(epsilon)
        partials[2,2] = (DataArray[_fz, iBlockP, i  , j  , 1] - DataArray[_fz, iBlockP, i  , j  , 0])/(epsilon)
    elif k == nK-1:
        partials[2,0] = (DataArray[_fx, iBlockP, i  , j  , nK-1] - DataArray[_fx, iBlockP, i  , j  , nK-2])/(epsilon)
        partials[2,1] = (DataArray[_fy, iBlockP, i  , j  , nK-1] - DataArray[_fy, iBlockP, i  , j  , nK-2])/(epsilon)
        partials[2,2] = (DataArray[_fz, iBlockP, i  , j  , nK-1] - DataArray[_fz, iBlockP, i  , j  , nK-2])/(epsilon)
    else:
        partials[2,0] = (DataArray[_fx, iBlockP, i  , j  , k+1] - DataArray[_fx, iBlockP, i  , j  , k-1])/(2*epsilon)
        partials[2,1] = (DataArray[_fy, iBlockP, i  , j  , k+1] - DataArray[_fy, iBlockP, i  , j  , k-1])/(2*epsilon)
        partials[2,2] = (DataArray[_fz, iBlockP, i  , j  , k+1] - DataArray[_fz, iBlockP, i  , j  , k-1])/(2*epsilon)

    return partials
