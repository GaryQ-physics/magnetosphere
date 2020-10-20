import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import util
import regions

run = 'DIPTSUR2'
cut = False
para = True

if os.path.exists('/home/gary/'):
    pointfile_path = '/home/gary/Downloads/points_for_gary-short.txt'
else:
    pointfile_path = '/home/gquaresi/points_for_gary.txt'


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

if not os.path.exists(direct):
    os.makedirs(direct)

pm = 31.875
reg =  {'xlims': (-pm, pm),
        'ylims': (-pm, pm),
        'zlims': (-pm, pm),
        'd': 0.25
        }

points = np.loadtxt(pointfile_path)
#results = np.nan*np.empty((points.shape[0],3))

def RUN(i):
    point = points[i,:]
    reg2 = reg.copy()

    def regionconfig(i, lims)
        if point[i] > 0:
            minn = reg['d']*round((point[i] - pm)/reg['d'])
            maxx = reg['d']*round((point[i] + pm)/reg['d'])
            minn = max(minn, pm + reg['d'])
            reg2[lims] = (minn, maxx)
        if point[i] < 0:
            minn = reg['d']*round((point[i] - pm)/reg['d'])
            maxx = reg['d']*round((point[i] + pm)/reg['d'])
            maxx = min(maxx, -pm - reg['d'])
            reg2[lims] = (minn, maxx)

    regionconfig(0, 'xlims')
    regionconfig(1, 'ylims')
    regionconfig(2, 'zlims')


    result = regions.signedintegrate(run, time, points[i,:], regions=(reg,reg2), rmin=rmin, locationtype='GSM')

    assert(len(result.shape)==3)
    if result.shape[0]==1:
        result = result[0, 2, :]
    else:
        result = np.sum(result, axis=0)


    print('i=%d_results_'%(i) + str(result[0]) + ',' + str(result[1]) + ',' + str(result[2]) + ',' + '\n')

    return result


if para:
    from joblib import Parallel, delayed
    import multiprocessing
    num_cores = multiprocessing.cpu_count()
    if num_cores is not None and num_cores > points.shape[0] :
        num_cores = points.shape[0]
    print('Parallel processing {0:d} point(s) using {1:d} cores'\
          .format(points.shape[0], num_cores))
    results = Parallel(n_jobs=num_cores)(\
            delayed(RUN)(i) for i in range(points.shape[0]))
else:
    results = []
    for i in range(points.shape[0]):
        results.append(RUN(i))

results = np.array(results)


bs_results_nopoints = open(direct + 'bs_results_nopoints.csv','w')
bs_results = open(direct + 'bs_results.csv','w')
for i in range(points.shape[0]):
    bs_results_nopoints.write(str(results[i,0]) + ',')
    bs_results_nopoints.write(str(results[i,1]) + ',')
    bs_results_nopoints.write(str(results[i,2]) + ',')

    bs_results.write(str(points[i,0]) + ',')
    bs_results.write(str(points[i,1]) + ',')
    bs_results.write(str(points[i,2]) + ',')
    bs_results.write(str(results[i,0]) + ',')
    bs_results.write(str(results[i,1]) + ',')
    bs_results.write(str(results[i,2]) + ',')
bs_results_nopoints.close()
bs_results.close()

txt = open(direct + 'bs_results.txt', 'w')
txt.write('time = (%d,%d,%d,%d,%d,%d,%d), rmin=%f, run=%s\n'%(time + (rmin,run)))
txt.write(str(reg)+'\n')
txt.write('point_x point_y point_z deltaB_x deltaB_x deltaB_x\n')
np.savetxt(txt, np.column_stack([points, results]), fmt='%.5f')
txt.close()

