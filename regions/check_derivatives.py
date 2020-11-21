import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
from probe import GetRunData
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

def GetFrobeniusNormDel(delF):
    ret = np.nan*np.empty(delF.shape[0])
    for i in range(delF.shape[0]):
        ret[i,:] = np.linalg.norm(delF[i,:,:], ord='fro') # frobenius norm, equivalently 2-norm of flattened vector. It's equivalent to set ord=None
    return ret

def GetOperatorNormDel(delF):
    ret = np.nan*np.empty(delF.shape[0])
    for i in range(delF.shape[0]):
        ret[i,:] = np.linalg.norm(delF[i,:,:], ord=2) # operator norm inherited from column vector 2-norm, equivalently largest singular value.
    return ret



####################
run = 'DIPTSUR2'
cut = True
pntlist = 'native_random_sampled'
####################

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
if cut:
    rmin = rCurrents
    direct = direct + 'excluding_currents_before_rCurrents/'
else:
    rmin = 0.
    direct = direct + 'including_currents_before_rCurrents/'

points = np.fromfile(direct + 'derivatives_points.bin').reshape((273,3))
results = np.fromfile(direct + 'derivatives_results.bin').reshape((3,273,3,3))

del_j_batsrus = results[0,:,:,:]
del_b_batsrus = results[1,:,:,:]
del_b1_batsrus = results[2,:,:,:]

div_Bbats = GetDivergence(del_b_batsrus)
div_B1bats = GetDivergence(del_b1_batsrus)
div_Jbats = GetDivergence(del_j_batsrus)

curl_Bbats = GetCurl(del_b_batsrus)
curl_B1bats = GetCurl(del_b1_batsrus)
curl_Jbats = GetCurl(del_j_batsrus)

FrobeniusNormDel_Bbats = GetFrobeniusNormDel(del_b_batsrus)
FrobeniusNormDel_B1bats = GetFrobeniusNormDel(del_b1_batsrus)
FrobeniusNormDel_Jbats = GetFrobeniusNormDel(del_j_batsrus)

OperatorNormDel_Bbats = GetOperatorNormDel(del_b_batsrus)
OperatorNormDel_B1bats = GetOperatorNormDel(del_b1_batsrus)
OperatorNormDel_Jbats = GetOperatorNormDel(del_j_batsrus)

Bbats = GetRunData(run, time, points, 'b')
B1bats = GetRunData(run, time, points, 'b1')
Jbats = GetRunData(run, time, points, 'j')

unitmu0 = phys['mu0'] * (phys['muA']/(phys['m']**2))
print('unitmu0 = %f'%(unitmu0))

'''
f = open(direct + 'gauss_check.txt','w')
f.write('point_x point_y point_z divB_sim divB1_sim divJ_sim\n')
np.savetxt(f, np.column_stack([points, GetDivergence(del_b_batsrus), GetDivergence(del_b1_batsrus), GetDivergence(del_j_batsrus)]))
f.close()


f = open(direct + 'ampere_check.txt','w')
f.write('point_x point_y point_z curlBx_sim curlBy_sim curlBz_sim Jx_scaled Jy_scaled Jz_scaled\n')
np.savetxt(f, np.column_stack([points, GetCurl(del_b_batsrus), J_scaled]))
f.close()
'''

outname = direct + 'derivatives_df.txt'
print('writing ' + outname)
f = open(outname, 'w')
f.write('point_x point_y point_z')
f.write(' div_Bbats div_B1bats div_Jbats')
f.write(' curl_Bbats_x curl_Bbats_y curl_Bbats_z')
f.write(' curl_B1bats_x curl_B1bats_y curl_B1bats_z')
f.write(' Bbats_x Bbats_y Bbats_z')
f.write(' B1bats_x B1bats_y B1bats_z')
f.write(' Jbats_x Jbats_y Jbats_z')
f.write(' FrobeniusNormDel_Bbats FrobeniusNormDel_B1bats FrobeniusNormDel_Jbats')
f.write(' OperatorNormDel_Bbats OperatorNormDel_B1bats OperatorNormDel_Jbats')
f.write('\n')

arr = np.column_stack([ div_Bbats, div_B1bats, div_Jbats,
                        curl_Bbats_x, curl_Bbats_y, curl_Bbats_z,
                        curl_B1bats_x, curl_B1bats_y, curl_B1bats_z,
                        Bbats_x, Bbats_y, Bbats_z,
                        B1bats_x, B1bats_y, B1bats_z,
                        Jbats_x, Jbats_y, Jbats_z,
                        FrobeniusNormDel_Bbats, FrobeniusNormDel_B1bats, FrobeniusNormDel_Jbats,
                        OperatorNormDel_Bbats, OperatorNormDel_B1bats, OperatorNormDel_Jbats ])

np.savetxt(f, arr)
f.close()
