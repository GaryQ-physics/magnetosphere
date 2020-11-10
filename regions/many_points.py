import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import util
import regions
import dissection as di

run = 'DIPTSUR2'
cut = True
para = True

if os.path.exists('/home/gary/'):
    pointfile_path = '/home/gary/Downloads/points_for_gary-test.txt'
else:
    pointfile_path = '/home/gquaresi/yplus_points.txt'


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

if not os.path.exists(direct):
    os.makedirs(direct)

points = np.loadtxt(pointfile_path)


#results = np.nan*np.empty((points.shape[0],3))

#pm = 31.875
#reg =  {'xlims': (-pm, pm),
#        'ylims': (-pm, pm),
#        'zlims': (-pm, pm),
#        'd': 0.25
#        }
def RUN(i):
    print('i = %d_start_'%(i))
    regs = di.GetRegions(points[i,:])
    result = regions.signedintegrate(run, time, points[i,:], regions=regs, rmin=rmin, locationtype='GSM')

    assert(len(result.shape)==3 and result.shape[1]==3 and result.shape[2]==3)
    if result.shape[0]==1:
        result = result[0, 2, :]
    else:
        result = np.sum(result, axis=0)[2,:]

    print('i = %d_results_'%(i) + str(result[0]) + ',' + str(result[1]) + ',' + str(result[2]) + ',' + '\n')

    return result

import time as tm
t0 = tm.time()

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


print('writing ' + direct + 'bs_results_nopoints.csv')
bs_results_nopoints = open(direct + 'bs_results_nopoints.csv','w')
print('writing ' + direct + 'bs_results.csv')
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
print('wrote ' + direct + 'bs_results_nopoints.csv')
bs_results.close()
print('wrote ' + direct + 'bs_results.csv')

print('writing ' + direct + 'bs_results.txt')

util.safeprep_fileout(direct + 'bs_results.txt')
txt = open(direct + 'bs_results.txt', 'w')

txt.write('time = (%d,%d,%d,%d,%d,%d,%d), rmin=%f, run=%s, cut=%s\n\n'%(time + (rmin,run,cut)))
txt.write('point_x point_y point_z Bx_biotsavart By_biotsavart Bz_biotsavart\n')
np.savetxt(txt, np.column_stack([points, results]), fmt='%.5f')
txt.close()
print('wrote ' + direct + 'bs_results.txt')

print('writing ' + direct + 'corresponding_regions.txt')
regstxt = open(direct + 'corresponding_regions.txt', 'w')
for i in range(points.shape[0]):
    regs = di.GetRegions(points[i,:])
    regstxt.write(str(regs))
    regstxt.write('\n\n')
regstxt.close()
print('wrote ' + direct + 'corresponding_regions.txt')



print('many_points ran in %f hours'%((tm.time()-t0)/3600.))#ran in 3.415081 hours
