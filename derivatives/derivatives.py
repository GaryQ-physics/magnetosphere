import os
import numpy as np

from config import conf
from probe import GetRunData
from units_and_constants import phys
import biot_savart as bs
import dissection as di


def GetDel(run, time, field, points, epsilon=0.0625, para=False, debug=False):
    # field = 'b_biotsavart' ; 'b_batsrus' ; 'b1_batsrus' ; 'j_batsrus' ; 'testinterp'

    def func(i):
        point = points[i,:].copy()
        xplus = point + epsilon*np.array([1,0,0])
        xmin  = point - epsilon*np.array([1,0,0])
        yplus = point + epsilon*np.array([0,1,0])
        ymin  = point - epsilon*np.array([0,1,0])
        zplus = point + epsilon*np.array([0,0,1])
        zmin  = point - epsilon*np.array([0,0,1])
        pts = np.array([xplus, xmin, yplus, ymin, zplus, zmin])
        if debug: print(pts.shape)

        if field == 'b_biotsavart':
            #pm = 31.875
            #reg =  {'xlims': (-pm, pm),
            #        'ylims': (-pm, pm),
            #        'zlims': (-pm, pm),
            #        'd': 0.25
            #        }
            regs = di.GetRegions(point)
            Fs = bs.biot_savart_run(run, time, pts, regs, summed=True, separateRegions=False)
        elif field == 'testinterp':
            Fs = np.column_stack([Bx_testinterp(pts), By_testinterp(pts), Bz_testinterp(pts)])
        elif '_batsrus' in field:
            fi = field.split('_batsrus')[0]
            if debug: print('fi=%s'%(fi))
            Fs = GetRunData(run, time, pts, fi)
        else:
            raise ValueError('invalid field string')

        #print(Fs.shape)
        delF = np.nan*np.empty((3,3))
        delF[0,:] = ( Fs[0,:] - Fs[1,:] )/(2.*epsilon)
        delF[1,:] = ( Fs[2,:] - Fs[3,:] )/(2.*epsilon)
        delF[2,:] = ( Fs[4,:] - Fs[5,:] )/(2.*epsilon)
        #print(delF)

        return delF

    if para:
        assert(field != 'b_biotsavart') #would run out of memory anyway

        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is not None and num_cores > points.shape[0]:
            num_cores = points.shape[0]
        if debug: print('Parallel processing {0:d} point(s) using {1:d} cores'\
              .format(points.shape[0], num_cores))
        ret = Parallel(n_jobs=num_cores)(\
                delayed(func)(i) for i in range(points.shape[0]))

    else:
        if debug: print('Serial processing {0:d} point(s)'\
              .format(points.shape[0],))
        ret = []
        for i in range(points.shape[0]):
            ret.append(func(i))

    return np.array(ret)


def GetDel_vectorized(run, time, field, points, epsilon=0.0625, para=False, debug=False,library='kameleon'):

    pts = np.empty((points.shape[0], 6, 3))
    for i in range(points.shape[0]):
        pts[i,0,:] = points[i,:] + epsilon*np.array([1,0,0])
        pts[i,1,:] = points[i,:] - epsilon*np.array([1,0,0])
        pts[i,2,:] = points[i,:] + epsilon*np.array([0,1,0])
        pts[i,3,:] = points[i,:] - epsilon*np.array([0,1,0])
        pts[i,4,:] = points[i,:] + epsilon*np.array([0,0,1])
        pts[i,5,:] = points[i,:] - epsilon*np.array([0,0,1])

    assert('_batsrus' in field)
    if para:
        print('para is irrelevantly True')
    else:
        print('para is irrelevantly False')


    fi = field.split('_batsrus')[0]
    if debug: print('fi=%s'%(fi))

    if debug: print(pts)
    pts = pts.reshape((6*points.shape[0], 3))
    Fs = GetRunData(run, time, pts, fi, library=library).reshape((points.shape[0], 6, 3))
    pts = pts.reshape((points.shape[0], 6, 3))
    if debug: print(pts)

    delF = np.nan*np.empty((points.shape[0], 3, 3))
    for i in range(points.shape[0]):
        delF[i,0,:] = ( Fs[i,0,:] - Fs[i,1,:] )/(2.*epsilon)
        delF[i,1,:] = ( Fs[i,2,:] - Fs[i,3,:] )/(2.*epsilon)
        delF[i,2,:] = ( Fs[i,4,:] - Fs[i,5,:] )/(2.*epsilon)

    return delF


#def GetDivergence(delF):
#    divF = np.nan*np.empty(delF.shape[0])
#    for i in range(delF.shape[0]):
#        divF[i] = delF[i,0,0] + delF[i,1,1] + delF[i,2,2]
#
#    return divF

def GetDivergence(delF):
    return delF[:,0,0] + delF[:,1,1] + delF[:,2,2]


def GetCurl(delF):
    curlF = np.nan*np.empty((delF.shape[0], 3))
    for i in range(delF.shape[0]):
        curlF_tens = delF[i,:,:] - delF[i,:,:].transpose()
        curlF[i,:] = np.array([ curlF_tens[1,2], curlF_tens[2,0], curlF_tens[0,1] ])

    return curlF

def GetFrobeniusNormDel(delF):
    ret = np.nan*np.empty(delF.shape[0])
    for i in range(delF.shape[0]):
        ret[i] = np.linalg.norm(delF[i,:,:], ord='fro') # frobenius norm, equivalently 2-norm of flattened vector. It's equivalent to set ord=None
    return ret

def GetOperatorNormDel(delF):
    ret = np.nan*np.empty(delF.shape[0])
    for i in range(delF.shape[0]):
        ret[i] = np.linalg.norm(delF[i,:,:], ord=2) # operator norm inherited from column vector 2-norm, equivalently largest singular value.
    return ret



