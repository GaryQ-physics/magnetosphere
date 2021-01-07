# for derivatives.md
import os
import numpy as np

from config import conf
from probe import GetRunData
import util
from units_and_constants import phys
import biot_savart as bs
import dissection as di
from datetime import datetime
from derivatives import GetDel
from check_derivatives import GetDivergence, GetCurl, GetFrobeniusNormDel, GetOperatorNormDel

now = datetime.now()
log = open('derivatives_script-'+now.strftime("%Y%m%dT%H%M%S")+'.log', 'w')
log.write('script began' + now.strftime("%Y-%m-%d T%H:%M:%S") + '\n')
log.write('current working directory      '  + os.getcwd() + '\n')

####################
run = 'DIPTSUR2'
cut = True
pntlist = 'native_random_sampled'
debug = False
####################

###### calculating partial derivatives ######
if run == 'DIPTSUR2':
    time = (2019,9,2,6,30,0,0)
    rCurrents = 1.8
    rBody = 1.5
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
log.write('loading points from ' + direct + pntlist + '_points.txt')

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

results = np.nan*np.empty((3, points.shape[0], 3, 3))

if debug: print('computing arrays')

types = ['j_batsrus', 'b_batsrus', 'b1_batsrus']
epsilon=0.0625
for i in range(3):
    results[i,:,:,:] = GetDel(run, time, types[i], points, epsilon=epsilon, debug=debug)
log.write('epsilon=%f\n'%(epsilon))

if debug: print('writing arrays')
results.tofile(direct + 'derivatives_results.bin')
points.tofile(direct + 'derivatives_points.bin')
log.write('wrote "results" array to ' + direct + 'derivatives_results.bin\n')
log.write('flattened from shape %s\n with first index indexing:\n'%(str(results.shape)))
[log.write('    %d  ->  %s\n'%(i,types[i])) for i in range(3)]

if debug: print('wrote arrays')
if debug: print(results)

log.write('derivatives computed in %f minutes\n'%((tm.time()-t0)/60.))
if debug: print('derivatives computed in %f minuts'%((tm.time()-t0)/60.))

###### ploting #############
import pandas as pd

imagedir = conf['base'] + 'images/' + run + '/'
if not os.path.exists(imagedir):
    os.mkdir(imagedir)

assert(types == ['j_batsrus', 'b_batsrus', 'b1_batsrus'])
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
log.write('unitmu0 = %f\n'%(unitmu0))
if debug: print('unitmu0 = %f'%(unitmu0))

outname = direct + 'derivatives_df.txt'
if debug: print('writing ' + outname)
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

arr = np.column_stack([ points,
                        div_Bbats, div_B1bats, div_Jbats,
                        curl_Bbats,
                        curl_B1bats,
                        Bbats,
                        B1bats,
                        Jbats,
                        FrobeniusNormDel_Bbats, FrobeniusNormDel_B1bats, FrobeniusNormDel_Jbats,
                        OperatorNormDel_Bbats, OperatorNormDel_B1bats, OperatorNormDel_Jbats ])

np.savetxt(f, arr)
f.close()
log.write('wrote all data in text file ' + outname + 'to be imported as pandas dataframe\n')
if debug: print('wrote ' + outname)

data = pd.read_csv(outname, sep=" ")
#import pickle
#with open('derpic.pkl', 'wb') as handle:
#    pickle.dump(data, handle)

###### plotting #####

    ###prep##
distance = np.sqrt(data['point_x']**2 + data['point_y']**2 + data['point_z']**2)
if debug:
    print('distance:')
    print(distance)

cut = distance > rBody
if debug:
    print('cut:')
    print(cut)

normB = np.sqrt(data['Bbats_x']**2 + data['Bbats_y']**2 + data['Bbats_z']**2)
normB1 = np.sqrt(data['B1bats_x']**2 + data['B1bats_y']**2 + data['B1bats_z']**2)
normJ = np.sqrt(data['Jbats_x']**2 + data['Jbats_y']**2 + data['Jbats_z']**2)

ON_div_Bbats = data['div_Bbats']/data['OperatorNormDel_Bbats']
ON_div_B1bats = data['div_B1bats']/data['OperatorNormDel_B1bats']
ON_div_Jbats = data['div_Jbats']/data['OperatorNormDel_Jbats']

FN_div_Bbats = data['div_Bbats']/data['FrobeniusNormDel_Bbats']
FN_div_B1bats = data['div_B1bats']/data['FrobeniusNormDel_B1bats']
FN_div_Jbats = data['div_Jbats']/data['FrobeniusNormDel_Jbats']

Re_div_Bbats = data['div_Bbats']/normB
Re_div_B1bats = data['div_B1bats']/normB1
Re_div_Jbats = data['div_Jbats']/normJ

    ###save png###

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

fig, axes = plt.subplots(figsize=(12,4), nrows=1, ncols=3, dpi=300)

axes[0].plot(distance[cut], data['div_Bbats'][cut], '.')
axes[0].xaxis.set_minor_locator(AutoMinorLocator())
axes[0].set_yscale('symlog', linthreshy=0.01)
axes[0].set_xlabel('distance from center [$R_E$]')
axes[0].set_ylabel('div(B) [$\\frac{nT}{R_E}$]')
axes[0].set_title('swmf units')

axes[1].plot(distance[cut], Re_div_Bbats[cut], '.')
axes[1].xaxis.set_minor_locator(AutoMinorLocator())
axes[1].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[1].set_xlabel('distance from center [$R_E$]')
axes[1].set_ylabel('$\\frac{div(B)}{norm(B)}$ [$\\frac{1}{R_E}$]', fontsize=12)
axes[1].set_title('relative divergence')

axes[2].plot(distance[cut], ON_div_Bbats[cut], '.')
axes[2].xaxis.set_minor_locator(AutoMinorLocator())
axes[2].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[2].set_xlabel('distance from center [$R_E$]')
axes[2].set_ylabel('$\\frac{div(B)}{NormD(B)}$', fontsize=12)
axes[2].set_title('normalize divergence')

fig.suptitle('Divergence of B', fontsize=16)
#https://stackoverflow.com/questions/8248467/matplotlib-tight-layout-doesnt-take-into-account-figure-suptitle/45161551#45161551
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

fig.savefig(imagedir+'divergence_B.png')
if debug: print('saved png ' + imagedir+'divergence_B.png')
log.write('saved png ' + imagedir+'divergence_B.png\n')

#plt.close()
del axes
del fig
#del matplotlib
#del plt
#del AutoMinorLocator
#import matplotlib
#matplotlib.use("Agg")
#import matplotlib.pyplot as plt
#from matplotlib.ticker import AutoMinorLocator

fig, axes = plt.subplots(figsize=(12,4), nrows=1, ncols=3, dpi=300)

axes[0].plot(distance[cut], data['div_Jbats'][cut], '.')
axes[0].xaxis.set_minor_locator(AutoMinorLocator())
axes[0].set_yscale('symlog', linthreshy=0.01)
axes[0].set_xlabel('distance from center [$R_E$]')
axes[0].set_ylabel('div(J) [$\\frac{\mu A / m^2}{R_E}$]')
axes[0].set_title('swmf units')

axes[1].plot(distance[cut], Re_div_Jbats[cut], '.')
axes[1].xaxis.set_minor_locator(AutoMinorLocator())
axes[1].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[1].set_xlabel('distance from center [$R_E$]')
axes[1].set_ylabel('$\\frac{div(J)}{norm(J)}$ [$\\frac{1}{R_E}$]', fontsize=12)
axes[1].set_title('relative divergence')

axes[2].plot(distance[cut], ON_div_Jbats[cut], '.')
axes[2].xaxis.set_minor_locator(AutoMinorLocator())
axes[2].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[2].set_xlabel('distance from center [$R_E$]')
axes[2].set_ylabel('$\\frac{div(J)}{NormD(J)}$', fontsize=12)
axes[2].set_title('normalize divergence')

fig.suptitle('Divergence of J', fontsize=16)
#https://stackoverflow.com/questions/8248467/matplotlib-tight-layout-doesnt-take-into-account-figure-suptitle/45161551#45161551
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

fig.savefig(imagedir+'divergence_J.png')
if debug: print('saved png ' + imagedir+'divergence_J.png')
log.write('saved png ' + imagedir+'divergence_J.png\n')


now = datetime.now()
log.write('script ended' + now.strftime("%Y-%m-%d T%H:%M:%S") + '\n')
log.close()
