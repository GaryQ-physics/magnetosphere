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

# B = [0.5*y**2, z*x, z]
# curlB = [-x, 0, z-y]
# divB = 1.

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

def divBtest(X):
    ret = np.zeros(X.shape[0])
    for i in range(X.shape[0]):
        ret[i] = 1.
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


def GetCurlB(points, filename, method='biotsavart', para=False):
    pm = 31.875
    reg =  {'xlims': (-pm, pm),
            'ylims': (-pm, pm),
            'zlims': (-pm, pm),
            'd': 0.25
            }
    epsilon = 1./16.

    def func(i):
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
            Bs = bs.biot_savart_run(run, time, pts, reg)
        if method=='b_batsrus':
            Bs = probe(filename, pts, var=['bx','by','bz'], library='kameleon')
        if method=='b1_batsrus':
            Bs = probe(filename, pts, var=['b1x','b1y','b1z'], library='kameleon')
        if method=='test':
            Bs = Btest(pts)
        if method=='testinterp':
            Bs = np.column_stack([Bx_testinterp(pts), By_testinterp(pts), Bz_testinterp(pts)])

        #print(Bs.shape)
        delB = np.nan*np.empty((3,3))
        delB[0,:] = ( Bs[0,:] - Bs[1,:] )/(2.*epsilon)
        delB[1,:] = ( Bs[2,:] - Bs[3,:] )/(2.*epsilon)
        delB[2,:] = ( Bs[4,:] - Bs[5,:] )/(2.*epsilon)
        #print(delB)

        curlB_tens = delB - delB.transpose()
        #print(curlB_tens)
        curlB = np.array([ curlB_tens[1,2], curlB_tens[2,0], curlB_tens[0,1] ])
        print(curlB)
        #print('\n\n\n')
        return curlB

    if para:
        assert(False) #would run out of memory anyway
    else:
        ret = []
        for i in range(points.shape[0]):
            ret.append(func(i))

    return np.array(ret)

run = 'DIPTSUR2'
method = 'b_batsrus'
cut = True
debug = False

if run == 'DIPTSUR2':
    time = (2019,9,2,6,30,0,0)
    rCurrents = 1.8
if run == 'IMP10_RUN_SAMPLE':
    time = (2019,9,2,7,0,0,0)
    rCurrents = 1.7
if run == 'TESTANALYTIC':
    time = (2000,1,1,0,10,0,0)
    rCurrents = 0.

direct = conf[run+'_derived'] + 'regions/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6)
if cut:
    rmin = rCurrents
    direct = direct + 'excluding_currents_before_rCurrents/'
else:
    rmin = 0.
    direct = direct + 'including_currents_before_rCurrents/'


if os.path.exists('/home/gary/'):
    #data = np.loadtxt(direct + 'including_currents_before_rCurrents/bs_results.txt', skiprows=3)
    #points = data[:, 0:3]

    #points = np.array([[5.2,4.3,-2.1],[2.,2.,7.]])
    #points = np.column_stack([np.arange(2,30),np.zeros(28),np.zeros(28)])

    points = np.loadtxt('/home/gary/Downloads/points_for_gary-test.txt')
else:
    points = np.loadtxt('/home/gquaresi/full_points.txt')

if debug:
    print(points)


filename = util.time2CDFfilename(run, time)

import time as tm
t0 = tm.time()


J = probe(filename, points, var=['jx','jy','jz'], library='kameleon')
J_scaled = phys['mu0'] * (phys['muA']/(phys['m']**2)) *J
curlB = GetCurlB(points, filename, method=method)

#error = np.einsum('ij,ij->i', curlB - J_scaled, curlB - J_scaled)
#error = error/np.einsum('ij,ij->i', J_scaled, J_scaled)
#error = np.sqrt(error)

if debug:
    print(J_scaled)
    print(curlB)
    print(curlBtest(points))
    #print(error)

if not os.path.exists(direct):
    os.makedirs(direct)

outfname = direct+'check_curl_method_%s.txt'%(method)
util.safeprep_fileout(outfname)
txt = open(outfname, 'w')

txt.write('time = (%d,%d,%d,%d,%d,%d,%d),  run=%s, (units: nT and R_e (J scaled by mu0))\n\n'%(time + (run,)))
txt.write('point_x point_y point_z Jx Jy Jz curlBx_%s curlBy_%s curlBz_%s\n'%(method,method,method))
np.savetxt(txt, np.column_stack([points, J_scaled, curlB]), fmt='%.5f')

txt.close()

print('runtime = %f minutes'%((tm.time()-t0)/60.)) #110 min biotsavart, 0.1 min batsrus

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

####### Oct 18th 2020 13:52 ########## method='biotsavart' #########
(python2.7) gquaresi@sunspot ~/magnetosphere/regions $ python check_curl.py 


 util path: ~/magnetosphere/util/util.py 


[[ 5.2  4.3 -2.1]
 [ 2.   2.   7. ]]
Kameleon::close() calling model's close
Kameleon::close() calling model's close
Kameleon::close() calling model's close
(6, 3)
Kameleon::close() calling model's close


hellothere


True
heh
(6, 3)
[[ 15.05304798  31.71843829   8.27645302]
 [ 28.81623777  12.54499775  -0.17297768]
 [  6.86216141   6.89162133 -27.59843898]]
[[ 0.          2.90220052  1.4142916 ]
 [-2.90220052  0.         -7.064599  ]
 [-1.4142916   7.064599    0.        ]]
[-7.064599   -1.4142916   2.90220052]




(6, 3)
Kameleon::close() calling model's close


hellothere


True
heh
(6, 3)
[[-71.32262607 -13.6614482   62.73757699]
 [-35.15499014  25.51301891  -1.07144681]
 [-87.70847478 -26.28448631  53.42347172]]
[[   0.           21.49354194  150.44605177]
 [ -21.49354194    0.           25.2130395 ]
 [-150.44605177  -25.2130395     0.        ]]
[  25.2130395  -150.44605177   21.49354194]




[[ -0.04388906   0.051706     0.22988535]
 [  5.34341367 -40.64481063  11.30561658]]
[[ 3.31596471e-02  2.07271695e-01 -1.14144844e+02]
 [ 1.61457733e+02  2.52562637e+01 -3.49682861e+02]]
[[ -71.14097595  -75.71589661 -165.93405151]
 [ 228.52600098   75.39697266 -247.07489014]]
[[  -7.064599     -1.4142916     2.90220052]
 [  25.2130395  -150.44605177   21.49354194]]
[[ -5.2    0.   -57.57]
 [ -2.     0.    -5.  ]]
runtime = 11.760381 minutes
#################

####### Oct 19th 2020 15:30 ########## method='biotsavart' #########
(python2.7) gary@gary-Inspiron-5567:~/magnetosphere/regions$ python check_curl.py 


 util path: ~/magnetosphere/util/util.py 


[[ 5.2  4.3 -2.1]
 [ 2.   2.   7. ]]
Kameleon::close() calling model's close
Kameleon::close() calling model's close
Kameleon::close() calling model's close
(6, 3)
Kameleon::close() calling model's close


 from deltaB 


True
checkpoint
(6, 3)
[[ 15.05304798  31.71843829   8.27645302]
 [ 28.81623777  12.54499775  -0.17297768]
 [  6.86216141   6.89162133 -27.59843898]]
[[ 0.          2.90220052  1.4142916 ]
 [-2.90220052  0.         -7.064599  ]
 [-1.4142916   7.064599    0.        ]]
[-7.064599   -1.4142916   2.90220052]




(6, 3)
Kameleon::close() calling model's close


 from deltaB 


True
checkpoint
(6, 3)
[[-71.32262607 -13.6614482   62.73757699]
 [-35.15499014  25.51301891  -1.07144681]
 [-87.70847478 -26.28448631  53.42347172]]
[[   0.           21.49354194  150.44605177]
 [ -21.49354194    0.           25.2130395 ]
 [-150.44605177  -25.2130395     0.        ]]
[  25.2130395  -150.44605177   21.49354194]




[[ -0.04388906   0.051706     0.22988535]
 [  5.34341367 -40.64481063  11.30561658]]
[[ 3.31596471e-02  2.07271695e-01 -1.14144844e+02]
 [ 1.61457733e+02  2.52562637e+01 -3.49682861e+02]]
[[ -71.14097595  -75.71589661 -165.93405151]
 [ 228.52600098   75.39697266 -247.07489014]]
[[  -7.064599     -1.4142916     2.90220052]
 [  25.2130395  -150.44605177   21.49354194]]
[[ -5.2    0.   -57.57]
 [ -2.     0.    -5.  ]]
[31.93329888  2.63489412]
runtime = 5.335643 minutes
#################
'''


