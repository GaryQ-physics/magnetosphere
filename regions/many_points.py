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
    pointfile_path = '/home/gary/Downloads/points_for_gary-test.txt'
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
    result = regions.signedintegrate(run, time, points[i,:], regions=(reg,), rmin=rmin, locationtype='GSM')

    assert(result.shape[0]==1 and len(result.shape)==3)
    result = result[0, 2, :]

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


csv_results = open(direct + 'csv_results.csv','w')
csv_full = open(direct + 'csv_full.csv','w')
for i in range(points.shape[0]):
    csv_results.write(str(results[i,0]) + ',')
    csv_results.write(str(results[i,1]) + ',')
    csv_results.write(str(results[i,2]) + ',')

    csv_full.write(str(points[i,0]) + ',')
    csv_full.write(str(points[i,1]) + ',')
    csv_full.write(str(points[i,2]) + ',')
    csv_full.write(str(results[i,0]) + ',')
    csv_full.write(str(results[i,1]) + ',')
    csv_full.write(str(results[i,2]) + ',')
csv_results.close()
csv_full.close()

txt_full = open(direct + 'txt_full.txt', 'w')
txt_full.write('time = (%d,%d,%d,%d,%d,%d,%d), rmin=%f, run=%s\n'%(time + (rmin,run)))
txt_full.write(str(reg)+'\n')
txt_full.write('point_x point_y point_z deltaB_x deltaB_x deltaB_x\n')
np.savetxt(txt_full, np.column_stack([points, results]), fmt='%.5f')
txt_full.close()

