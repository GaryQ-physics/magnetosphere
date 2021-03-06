import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
from probe import probe
import util
from units_and_constants import phys
import biot_savart as bs
import dissection as di

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



def GetDel(run, time, field, points, para=False, epsilon=0.0625, debug=False):
    filename = util.time2CDFfilename(run, time)

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
            Fs = probe(filename, pts, var=[fi+'x',fi+'y',fi+'z'], library='kameleon')
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
        assert(False) #would run out of memory anyway
    else:
        ret = []
        for i in range(points.shape[0]):
            ret.append(func(i))

    return np.array(ret)

if __name__=='__main__':
    print('RUNNING regions/derivates.py')

    ####################
    run = 'DIPTSUR2'
    cut = True
    pntlist = 'native_random_sampled'
    ####################
    debug = False

    if run == 'DIPTSUR2':
        time = (2019,9,2,6,30,0,0)
        rCurrents = 1.8
    if run == 'IMP10_RUN_SAMPLE':
        time = (2019,9,2,7,0,0,0)
        rCurrents = 1.7
    if run == 'TESTANALYTIC':
        time = (2000,1,1,0,10,0,0)
        rCurrents = 1.5

    direct = conf[run+'_derived'] + 'regions/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6)
    direct = direct + pntlist + '/'

    points = np.loadtxt(direct + pntlist + '_points.txt')
    if debug: print(points)

    if cut:
        rmin = rCurrents
        direct = direct + 'excluding_currents_before_rCurrents/'
    else:
        rmin = 0.
        direct = direct + 'including_currents_before_rCurrents/'

    if not os.path.exists(direct):
        os.mkdir(direct)

    import time as tm
    t0 = tm.time()

    if False:
        results = np.nan*np.empty((3, points.shape[0], 3, 3))

        results[0,:,:,:] = GetDel(run, time, 'j_batsrus', points)
        results[1,:,:,:] = GetDel(run, time, 'b_batsrus', points)
        results[2,:,:,:] = GetDel(run, time, 'b1_batsrus', points)

        print('writing arrays')
        results.tofile(direct + 'derivatives_results.bin')
        points.tofile(direct + 'derivatives_points.bin')

        if debug: print(results)
    else:
        delBbs = GetDel(run, time, 'b_biotsavart', points)
        print('writing array')
        delBbs.tofile(direct + 'derivatives_bs.bin')

    print('derivatives ran in %f minuts'%((tm.time()-t0)/60.))
