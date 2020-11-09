import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
from probe import probe
import util
from units_and_constants import phys

def GetDivergence(delF):
    divF = np.nan*np.empty(delF.shape[0])
    for i in range(delF.shape[0]):
        divF[i] = delF[i,0,0] + delF[i,1,1] + delF[i,2,2]

    return divF

def GetCurl(delF):
    curlF = np.nan*np.empty((delF.shape[0], 3))
    for i in range(delF.shape[0]):
        curlF_tens = delF[i,:,:] - delF[i,:,:].transpose()
        curlF[i,:] = np.array([ curlF_tens[1,2], curlF_tens[2,0], curlF_tens[0,1] ])

    return curlF

run = 'DIPTSUR2'
cut = True

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

filename = util.time2CDFfilename(run, time)


points = np.fromfile(direct + 'derivatives_points.bin').reshape((273,3))
results = np.fromfile(direct + 'derivatives_results.bin').reshape((3,273,3,3))

del_j_batsrus = results[0,:,:,:]
del_b_batsrus = results[1,:,:,:]
del_b1_batsrus = results[2,:,:,:]


f = open(direct + 'gauss_check.txt','w')
f.write('point_x point_y point_z divB_sim divB1_sim, divJ_sim')
np.savetxt(f, np.column_stack([points, GetDivergence(del_b_batsrus), GetDivergence(del_b1_batsrus), GetDivergence(del_j_batsrus)]))
f.close()

J = probe(filename, points, var=['jx','jy','jz'], library='kameleon')
J_scaled = phys['mu0'] * (phys['muA']/(phys['m']**2)) * J

f = open(direct + 'ampere_check.txt','w')
f.write('point_x point_y point_z curlBx_sim curlBy_sim curlBz_sim Jx_scaled Jy_scaled Jz_scaled')
np.savetxt(f, np.column_stack([points, GetCurl(del_b_batsrus), J_scaled]))
f.close()
