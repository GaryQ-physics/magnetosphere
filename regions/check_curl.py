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

# B = [0.5*y**2, z*x, z]
# curlB = [-x, 0, z-y]

def Btest(X):
    ret = np.zeros(X.shape)
    for i in range(X.shape[0]):
        #ret[i,:] = np.array([-X[i,1], X[i,0], 0.])
        ret[i,:] = np.array([X[i,1]**3, X[i,2]*X[i,0], X[i,2]]) #central difference quotient is exact for quadratics
    return ret

def curlBtest(X):
    ret = np.zeros(X.shape)
    for i in range(X.shape[0]):
        #ret[i,:] = np.array([0,0,2])
        ret[i,:] = np.array([-X[i,0], 0., X[i,2]-3.*X[i,1]**2])
    return ret

#adapted from magnetosphere/cutplane/fieldlines.py
from scipy.interpolate import RegularGridInterpolator
xax = np.linspace(-30, 30, 61)
yax = np.linspace(-30, 30, 61)
zax = np.linspace(-30, 30, 61)
Nx = xax.size
Ny = yax.size
Nz = zax.size

G2, G1, G3 = np.meshgrid(yax, xax, zax) # different than in make_grid
P = np.column_stack( (G1.flatten(), G2.flatten(), G3.flatten()) )

res = Btest(P)
Bx_interdata = res[:,0].reshape(Nx, Ny, Nz)
By_interdata = res[:,1].reshape(Nx, Ny, Nz)
Bz_interdata = res[:,2].reshape(Nx, Ny, Nz)

Bx_testinterp = RegularGridInterpolator((xax,yax,zax), Bx_interdata)
By_testinterp = RegularGridInterpolator((xax,yax,zax), By_interdata)
Bz_testinterp = RegularGridInterpolator((xax,yax,zax), Bz_interdata)


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
        point = points[i,:].copy()
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
        if method=='b_batsrus':
            Bs = probe(filename, pts, var=['bx','by','bz'], library='kameleon')
        if method=='b1_batsrus':
            Bs = probe(filename, pts, var=['b1x','b1y','b1z'], library='kameleon')
        if method=='test':
            Bs = Btest(pts)
        if method=='testinterp':
            Bs = np.column_stack([Bx_testinterp(pts), By_testinterp(pts), Bz_testinterp(pts)])

        print(Bs.shape)
        delB = np.nan*np.empty((3,3))
        delB[0,:] = ( Bs[0,:] - Bs[1,:] )/(2.*epsilon)
        delB[1,:] = ( Bs[2,:] - Bs[3,:] )/(2.*epsilon)
        delB[2,:] = ( Bs[4,:] - Bs[5,:] )/(2.*epsilon)
        print(delB)

        curlB_tens = delB - delB.transpose()
        print(curlB_tens)
        curlB = np.array([ curlB_tens[1,2], curlB_tens[2,0], curlB_tens[0,1] ])
        print(curlB)

        ret.append(curlB)

        print('\n\n\n')
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
points = np.array([[5.2,4.3,-2.1],[2.,2.,7.]])

print(points)


filename = util.time2CDFfilename(run, time)

import time as tm
t0 = tm.time()

B1 = probe(filename, points, var=['b1x','b1y','b1z'], library='kameleon')
B = probe(filename, points, var=['bx','by','bz'], library='kameleon')
J = probe(filename, points, var=['jx','jy','jz'], library='kameleon')*(phys['muA']/(phys['m']**2))
J_scaled = phys['mu0'] * J
curlB = GetCurlB(points, filename, method='testinterp')
print(J_scaled)
print(B)
print(B1)
print(curlB)
print(curlBtest(points))

print('runtime = %f minutes'%((tm.time()-t0)/60.))

'''
method='biotsavart'
[[-0.00111117  0.07677685 -0.08971719]
 [-0.00454559  0.06671563 -0.07704238]]
[[ 9.66720283e-02  4.36473012e-01 -1.13908623e+02]
 [ 1.04701966e-01  4.45897877e-01 -1.13904495e+02]]
[[  -5.92487621    0.97990096 -122.53599548]
 [  -4.15499353    0.87281436 -123.13650513]]
[[ 0.1731481   0.59475596 -0.15250371]
 [ 0.16899027  0.57523569 -0.1242519 ]]
'''


