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

if os.path.exists('/home/gary/'):
    pointfile_path = '/home/gary/Downloads/points_for_gary.txt'
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


pm = 31.875
reg =  {'xlims': (-pm, pm),
        'ylims': (-pm, pm),
        'zlims': (-pm, pm),
        'd': 0.25
        }

points = np.loadtxt(pointfile_path)
results = np.nan*np.empty((points.shape[0],3))

if not os.path.exists(direct):
    os.makedirs(direct)
csv_results = open(direct + 'csv_results.csv','w')
csv_full = open(direct + 'csv_full.csv','w')

for i in range(points.shape[0]):
    result = regions.signedintegrate(run, time, points[i,:], regions=(reg,), rmin=rmin, locationtype='GSM')

    assert(result.shape[0]==1 and len(result.shape)==3)
    result = result[0, 2, :]

    results[i, :] = result

    csv_results.write(str(result[0]) + ',')
    csv_results.write(str(result[1]) + ',')
    csv_results.write(str(result[2]) + ',')

    csv_full.write(str(points[i,0]) + ',')
    csv_full.write(str(points[i,1]) + ',')
    csv_full.write(str(points[i,2]) + ',')
    csv_full.write(str(result[0]) + ',')
    csv_full.write(str(result[1]) + ',')
    csv_full.write(str(result[2]) + ',')

csv_results.close()
csv_full.close()

txt_full = open(direct + 'txt_full.txt', 'w')
txt_full.write('time = (%d,%d,%d,%d,%d,%d,%d), rmin=%f, run=%s\n'%(time + (rmin,run)))
txt_full.write(str(reg)+'\n')
txt_full.write('point_x point_y point_z deltaB_x deltaB_x deltaB_x\n')
np.savetxt(txt_full, np.column_stack([points, results]), fmt='%.5f')
txt_full.close()

