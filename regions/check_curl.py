import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
from probe import probe
import util
from units_and_constants import phys
#import biot_savart_kameleon_interpolated_grid as bsk
import biot_savart as bs
from make_grid import make_grid, make_axes





def bsint(filename, pts, region):
    ax_list = make_axes(region['xlims'], region['ylims'], region['zlims'], region['d'])
    G = make_grid(ax_list, slices=False)
    J = probe(filename, G, var=['jx','jy','jz'], library='kameleon')*(phys['muA']/(phys['m']**2))
    deltaBs = bs.deltaB('deltaB', pts, G, J, V_char = region['d']**3)
    return deltaBs


def GetCurlB(points, filename, method='biotsavart'):
    pm = 31.875
    reg =  {'xlims': (-pm, pm),
            'ylims': (-pm, pm),
            'zlims': (-pm, pm),
            'd': 0.25
            }
    epsilon = 1./16.

    ret = []
    for i in range(points.shape[0]):
        point = points[i,:]
        xplus = point + epsilon*np.array([1,0,0])
        xmin =  point - epsilon*np.array([1,0,0])
        yplus = point + epsilon*np.array([0,1,0])
        ymin =  point - epsilon*np.array([0,1,0])
        zplus = point + epsilon*np.array([0,0,1])
        zmin =  point - epsilon*np.array([0,0,1])
        pts = np.array([xplus, xmin, yplus, ymin, zplus, zmin])
        print(pts.shape)

        if method=='biotsavart':
            Bs = bsint(filename, pts, reg)
        #Bs = probe(filename, pts, var=['bx','by','bz'], library='kameleon')
        #Bs = probe(filename, pts, var=['b1x','b1y','b1z'], library='kameleon')

        print(Bs.shape)
        delB = np.nan*np.empty((3,3))
        delB[0,:] = ( Bs[0,:] - Bs[1,:] )/(2.*epsilon)
        delB[1,:] = ( Bs[2,:] - Bs[3,:] )/(2.*epsilon)
        delB[2,:] = ( Bs[4,:] - Bs[5,:] )/(2.*epsilon)

        curlB_tens = delB - delB.transpose()
        curlB = np.array([ curlB_tens[1,2], curlB_tens[2,0], curlB_tens[0,1] ])
        ret.append(curlB)
    return np.array(ret)

run = 'DIPTSUR2'
cut = False

if run == 'DIPTSUR2':
    time = (2019,9,2,6,30,0,0)
    rCurrents = 1.8
if run == 'IMP10_RUN_SAMPLE':
    time = (2019,9,2,7,0,0,0)
    rCurrents = 1.7

direct = conf[run+'_derived'] + 'regions/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6)
if cut:
    rmin = rCurrents
    direct = direct + 'excluding_currents_before_rCurrents/'
else:
    rmin = 0.
    direct = direct + 'including_currents_before_rCurrents/'
data = np.loadtxt(direct + 'txt_full-short.txt', skiprows=3)
points = data[:, 0:3]
print(points)


filename = util.time2CDFfilename(run, time)

import time as tm
t0 = tm.time()

B1 = probe(filename, points, var=['b1x','b1y','b1z'], library='kameleon')
B = probe(filename, points, var=['bx','by','bz'], library='kameleon')
J = probe(filename, points, var=['jx','jy','jz'], library='kameleon')*(phys['muA']/(phys['m']**2))
J_scaled = phys['mu0'] * J
curlB = GetCurlB(points, filename, method='biotsavart')
print(J_scaled)
print(B)
print(B1)
print(curlB)
print('runtime = %f minutes'%((tm.time()-t0)/60.))
