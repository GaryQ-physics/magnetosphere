import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
from probe import probe
import util
from units_and_constants import phys

run = 'DIPTSUR2'
cut = True

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

filename = util.time2CDFfilename(run, time)

data = np.loadtxt(direct + 'bs_results.txt', skiprows=3)
points = data[:, 0:3]
dB = data[:, 3:6]
B1 = probe(filename, points, var=['b1x','b1y','b1z'], library='kameleon')
B = probe(filename, points, var=['bx','by','bz'], library='kameleon')

util.safeprep_fileout(direct + 'comparison.txt')
txt = open(direct + 'comparison.txt','w')

txt.write('B_bs is the field computed by biot savart integral of "j" in datafile\n')
txt.write('B_sim is the field given by "b" in datafile\n')
txt.write('B1_sim is the field given by "b1" in datafile\n\n')
txt.write('point_x point_y point_z Bx_bs By_bs Bz_bs B1x_sim B1y_sim B1z_sim Bx_sim By_sim Bz_sim\n')
np.savetxt(txt, np.column_stack([data, B1, B]), fmt='%.5f')
txt.close()
