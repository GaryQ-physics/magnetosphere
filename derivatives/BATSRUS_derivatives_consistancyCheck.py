import os
import sys
import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
from probe import GetRunData
import util
from units_and_constants import phys

from datetime import datetime
from derivatives import GetDel_vectorized, GetDivergence, GetCurl, GetFrobeniusNormDel, GetOperatorNormDel

now = datetime.now()
log = open('BATSRUS_derivatives_consistancyCheck-'+now.strftime("%Y%m%dT%H%M%S")+'.log', 'w')
log.write('script began' + now.strftime("%Y-%m-%d T%H:%M:%S") + '\n')
log.write('current working directory      '  + os.getcwd() + '\n')
log.write('USING VTK\n')
####################
run = 'LUHMANN1979'
pntlist = 'native_random_sampled2'
#pntlist = 'xz_plane_y=0.062500'
skip_computing = False
para = True
debug = False

library = 'vtk' #temporarily
log.write('using '+ library+' library')
####################

if run == 'DIPTSUR2':
    time = (2019,9,2,6,30,0,0)
    #time = (2019,9,2,4,10,0,0)
    rCurrents = 1.8
    rBody = 1.5

    epsilons = [1./16., 1./8.]
    corresponding_max_radii = [5., np.inf]

if run == 'IMP10_RUN_SAMPLE':
    time = (2019,9,2,7,0,0,0)
    rCurrents = 1.7
if run == 'TESTANALYTIC':
    time = (2000,1,1,0,10,0,0)
    rCurrents = 1.5
if run == 'LUHMANN1979':
    time = (2000,1,1,0,0,0,0)
    rCurrents = 1.1
    rBody = 1.

    epsilons = [1./16., 1./8.]
    corresponding_max_radii = [5., np.inf]

direct = conf[run+'_derived'] + '%s_library/'%(library) \
        + 'derivatives/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6) \
        + pntlist + '/'
if not os.path.exists(direct):
    os.makedirs(direct)

points_fname = conf['storage'] + 'pointlists/' + pntlist + '__points.txt'
points = np.loadtxt(points_fname)
if debug: print(points)

log.write('loaded points from ' + points_fname + '\n')
log.write('to plot with epsilons = %s and corresponding_max_radii = %s\n'%(str(epsilons),str(corresponding_max_radii)))

import time as tm

def GetDerivativesArray(epsil):
    results_fname = direct + 'partial_derivatives_epsilon=%f.bin'%(epsil)
    if not skip_computing:
        t0 = tm.time()
        results = np.nan*np.empty((3, points.shape[0], 3, 3))

        if debug: print('computing arrays')

        types = ['j_batsrus', 'b_batsrus', 'b1_batsrus']
        for i in range(3):
            results[i,:,:,:] = GetDel_vectorized(run, time, types[i], points, epsilon=epsil, para=para, debug=debug, library=library)
        log.write('epsilon=%f\n'%(epsil))

        if debug: print('writing arrays')
        results.tofile(results_fname)
        points.tofile(direct + 'derivatives_points.bin')
        if debug: print('wrote "results" array to ' + results_fname + '\n')
        log.write('wrote "results" array to ' + results_fname + '\n')
        log.write('flattened from shape %s\n with first index indexing:\n'%(str(results.shape)))
        [log.write('    %d  ->  %s\n'%(i,types[i])) for i in range(3)]

        if debug: print('wrote arrays')
        if debug: print(results)

        log.write('derivatives computed in %f minutes\n'%((tm.time()-t0)/60.))
        if debug: print('derivatives computed in %f minuts'%((tm.time()-t0)/60.))
    else:
        results = np.fromfile(results_fname).reshape((3, points.shape[0], 3, 3))
        if debug: print('loaded "results" array from ' + results_fname + '\n')
        log.write('loaded "results" array from ' + results_fname + '\n')

        pointsloaded = np.fromfile(direct + 'derivatives_points.bin').reshape((points.shape[0], 3))
        assert(np.all(points == pointsloaded))

    return results


###### process data #############
dels = [GetDerivativesArray(epsilon) for epsilon in epsilons]

DISTANCE = np.sqrt(points[:,0]**2 + points[:,1]**2 + points[:,2]**2)

imagedir = conf['base'] + 'images/' + run+'/' + 'using_%s/'%(library) \
        + 'derivatives/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6) \
        + pntlist + '/'
if not os.path.exists(imagedir):
    os.makedirs(imagedir)

#assert(types == ['j_batsrus', 'b_batsrus', 'b1_batsrus'])
del_j_batsrus = dels[0][0,:,:,:]
del_b_batsrus = dels[0][1,:,:,:]
del_b1_batsrus = dels[0][2,:,:,:]

for i in range(1, len(epsilons)):
    distanceSlice = np.logical_and(corresponding_max_radii[i-1] < DISTANCE,
                                    DISTANCE <= corresponding_max_radii[i])
    del_j_batsrus[distanceSlice,:,:] = dels[i][0,distanceSlice,:,:]
    del_b_batsrus[distanceSlice,:,:] = dels[i][1,distanceSlice,:,:]
    del_b1_batsrus[distanceSlice,:,:] = dels[i][2,distanceSlice,:,:]


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

Bbats = GetRunData(run, time, points, 'b', library=library)
B1bats = GetRunData(run, time, points, 'b1', library=library)
Jbats = GetRunData(run, time, points, 'j', library=library)

unitmu0 = phys['mu0'] * (phys['muA']/(phys['m']**2))
log.write('unitmu0 = %f\n'%(unitmu0))
if debug: print('unitmu0 = %f'%(unitmu0))

header = ''
header = header + 'point_x point_y point_z '
header = header + 'div_Bbats div_B1bats div_Jbats '
header = header + 'curl_Bbats_x curl_Bbats_y curl_Bbats_z '
header = header + 'curl_B1bats_x curl_B1bats_y curl_B1bats_z '
header = header + 'Bbats_x Bbats_y Bbats_z '
header = header + 'B1bats_x B1bats_y B1bats_z '
header = header + 'Jbats_x Jbats_y Jbats_z '
header = header + 'FrobeniusNormDel_Bbats FrobeniusNormDel_B1bats FrobeniusNormDel_Jbats '
header = header + 'OperatorNormDel_Bbats OperatorNormDel_B1bats OperatorNormDel_Jbats'

arr = np.column_stack([ points,
                        div_Bbats, div_B1bats, div_Jbats,
                        curl_Bbats,
                        curl_B1bats,
                        Bbats,
                        B1bats,
                        Jbats,
                        FrobeniusNormDel_Bbats, FrobeniusNormDel_B1bats, FrobeniusNormDel_Jbats,
                        OperatorNormDel_Bbats, OperatorNormDel_B1bats, OperatorNormDel_Jbats ])

data = pd.DataFrame(data=arr, columns=header.split(' '))

if True:
    data_fname = direct + 'derivatives_df.txt'
    if debug: print('writing ' + data_fname)
    f = open(data_fname, 'w')
    f.write(header)
    f.write('\n')
    np.savetxt(f, arr)
    f.close()
    log.write('wrote all data in text file ' + data_fname + ' to be imported as pandas dataframe\n')
    if debug: print('wrote ' + data_fname)

    dataload = pd.read_csv(data_fname, sep=" ")
    log.write('max diff columns:\n%s\n'%(str(np.max(np.abs(dataload-data)))))
    if np.all( np.abs(np.array(dataload-data)) > 1e-6 ):
        log.write('ERROR: abs > 1e-6')#!!!!!!!!


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

axes[0].plot(distance[cut], data['div_Jbats'][cut], '.')
#axes[0].xaxis.set_minor_locator(AutoMinorLocator())
axes[0].set_xscale('symlog', linthreshx=10., subsx=[2, 3, 4, 5, 6, 7, 8, 9])
axes[0].set_yscale('symlog', linthreshy=0.01)
axes[0].set_xlabel('distance from center [$R_E$]')
axes[0].set_ylabel('div(J) [$\\frac{\mu A / m^2}{R_E}$]')
axes[0].set_title('swmf units')

axes[1].plot(distance[cut], Re_div_Jbats[cut], '.')
#axes[1].xaxis.set_minor_locator(AutoMinorLocator())
axes[1].set_xscale('symlog', linthreshx=10., subsx=[2, 3, 4, 5, 6, 7, 8, 9])
axes[1].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[1].set_xlabel('distance from center [$R_E$]')
axes[1].set_ylabel('$\\frac{div(J)}{norm(J)}$ [$\\frac{1}{R_E}$]', fontsize=12)
axes[1].set_title('relative divergence')

axes[2].plot(distance[cut], ON_div_Jbats[cut], '.')
#axes[2].xaxis.set_minor_locator(AutoMinorLocator())
axes[2].set_xscale('symlog', linthreshx=10., subsx=[2, 3, 4, 5, 6, 7, 8, 9])
axes[2].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[2].set_xlabel('distance from center [$R_E$]')
axes[2].set_ylabel('$\\frac{div(J)}{NormD(J)}$', fontsize=12)
axes[2].set_title('normalize divergence')

fig.suptitle('Divergence of J', fontsize=16)
#https://stackoverflow.com/questions/8248467/matplotlib-tight-layout-doesnt-take-into-account-figure-suptitle/45161551#45161551
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

figname = 'divergence_J.png'
fig.savefig(imagedir+figname)
if debug: print('saved png ' + imagedir+figname)
log.write('saved png ' + imagedir+figname+'\n')

del axes
del fig

fig, axes = plt.subplots(figsize=(12,4), nrows=1, ncols=3, dpi=300)

axes[0].plot(distance[cut], data['div_Bbats'][cut], '.')
#axes[0].xaxis.set_minor_locator(AutoMinorLocator())
axes[0].set_xscale('symlog', linthreshx=10., subsx=[2, 3, 4, 5, 6, 7, 8, 9])
axes[0].set_yscale('symlog', linthreshy=0.01)
axes[0].set_xlabel('distance from center [$R_E$]')
axes[0].set_ylabel('div(B) [$\\frac{nT}{R_E}$]')
axes[0].set_title('swmf units')

axes[1].plot(distance[cut], Re_div_Bbats[cut], '.')
#axes[1].xaxis.set_minor_locator(AutoMinorLocator())
axes[1].set_xscale('symlog', linthreshx=10., subsx=[2, 3, 4, 5, 6, 7, 8, 9])
axes[1].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[1].set_xlabel('distance from center [$R_E$]')
axes[1].set_ylabel('$\\frac{div(B)}{norm(B)}$ [$\\frac{1}{R_E}$]', fontsize=12)
axes[1].set_title('relative divergence')

axes[2].plot(distance[cut], ON_div_Bbats[cut], '.')
#axes[2].xaxis.set_minor_locator(AutoMinorLocator())
axes[2].set_xscale('symlog', linthreshx=10., subsx=[2, 3, 4, 5, 6, 7, 8, 9])
axes[2].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[2].set_xlabel('distance from center [$R_E$]')
axes[2].set_ylabel('$\\frac{div(B)}{NormD(B)}$', fontsize=12)
axes[2].set_title('normalize divergence')

fig.suptitle('Divergence of B', fontsize=16)
#https://stackoverflow.com/questions/8248467/matplotlib-tight-layout-doesnt-take-into-account-figure-suptitle/45161551#45161551
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

figname = 'divergence_B.png'
fig.savefig(imagedir+figname)
if debug: print('saved png ' + imagedir+figname)
log.write('saved png ' + imagedir+figname+'\n')

del axes
del fig

fig, axes = plt.subplots(figsize=(12,4), nrows=1, ncols=3, dpi=300)

axes[0].plot(distance[cut], data['div_B1bats'][cut], '.')
#axes[0].xaxis.set_minor_locator(AutoMinorLocator())
axes[0].set_xscale('symlog', linthreshx=10., subsx=[2, 3, 4, 5, 6, 7, 8, 9])
axes[0].set_yscale('symlog', linthreshy=0.01)
axes[0].set_xlabel('distance from center [$R_E$]')
axes[0].set_ylabel('div(B1) [$\\frac{nT}{R_E}$]')
axes[0].set_title('swmf units')

axes[1].plot(distance[cut], Re_div_B1bats[cut], '.')
#axes[1].xaxis.set_minor_locator(AutoMinorLocator())
axes[1].set_xscale('symlog', linthreshx=10., subsx=[2, 3, 4, 5, 6, 7, 8, 9])
axes[1].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[1].set_xlabel('distance from center [$R_E$]')
axes[1].set_ylabel('$\\frac{div(B1)}{norm(B1)}$ [$\\frac{1}{R_E}$]', fontsize=12)
axes[1].set_title('relative divergence')

axes[2].plot(distance[cut], ON_div_B1bats[cut], '.')
#axes[2].xaxis.set_minor_locator(AutoMinorLocator())
axes[2].set_xscale('symlog', linthreshx=10., subsx=[2, 3, 4, 5, 6, 7, 8, 9])
axes[2].set_yscale('symlog', linthreshy=0.01, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
axes[2].set_xlabel('distance from center [$R_E$]')
axes[2].set_ylabel('$\\frac{div(B1)}{NormD(B1)}$', fontsize=12)
axes[2].set_title('normalize divergence')

fig.suptitle('Divergence of B1', fontsize=16)
#https://stackoverflow.com/questions/8248467/matplotlib-tight-layout-doesnt-take-into-account-figure-suptitle/45161551#45161551
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

figname = 'divergence_B1.png'
fig.savefig(imagedir+figname)
if debug: print('saved png ' + imagedir+figname)
log.write('saved png ' + imagedir+figname+'\n')

del axes
del fig

#### ampere using B ####
for cp in ['x', 'y', 'z']:
    fig, axes = plt.subplots(figsize=(12,4), nrows=1, ncols=2, dpi=300)

    axes[0].plot(distance, data['curl_Bbats_'+cp], '+', color='DarkRed', label='$curl(B)$')
    axes[0].plot(distance, unitmu0*data['Jbats_'+cp], 'x', color='Orange', label='$\\mu_0 J$')
    axes[0].set_yscale('symlog', linthreshy=1e-2)
    axes[0].legend()
    axes[0].set_title('dimensionfull values')
    axes[0].set_xlabel('distance from center [$R_E$]')
    axes[0].set_ylabel('values in  $\\frac{nT}{R_E}$')

    axes[1].plot(distance, (data['curl_Bbats_'+cp] - unitmu0*data['Jbats_'+cp])/(unitmu0*normJ),'x',
                 label='$\\frac{curl(B) - \mu_0 J}{\mu_0 norm(J)}$', color='LightBlue')
    axes[1].plot(distance, (data['curl_Bbats_'+cp] - unitmu0*data['Jbats_'+cp])/(unitmu0*data['Jbats_'+cp]),'+',
                 label='$\\frac{curl(B) - \mu_0 J}{\mu_0 J}$', color='Orange')
    axes[1].set_yscale('symlog', linthreshy=1e-4)
    axes[1].legend()
    axes[1].set_title('fractional errors')
    axes[1].set_xlabel('distance from center [$R_E$]')

    fig.suptitle(cp+' component', fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    figname = 'curlB_and_J_'+cp+'_component.png'
    fig.savefig(imagedir+figname)
    if debug: print('saved png ' + imagedir+figname)
    log.write('saved png ' + imagedir+figname+'n')

    del axes
    del fig

amp_err = np.sqrt( (data['curl_Bbats_x'] - unitmu0*data['Jbats_x'])**2
                  +(data['curl_Bbats_y'] - unitmu0*data['Jbats_y'])**2
                  +(data['curl_Bbats_z'] - unitmu0*data['Jbats_z'])**2 )

amp_percent_err = 100.*amp_err/(unitmu0*normJ)

fig, ax = plt.subplots(figsize=(12,4), dpi=300)

ax.plot(distance, amp_percent_err, '.')
ax.set_yscale('symlog', linthreshy=1e-2)
ax.vlines(rCurrents, 0, 100, linewidths=1.)
ax.set_title('percent error in amperes law')
ax.set_xlabel('distance from center [$R_E$]')
ax.set_ylabel('percent error (%)')

figname = 'curlB_and_J_percent_error.png'
fig.savefig(imagedir+figname)
if debug: print('saved png ' + imagedir+figname)
log.write('saved png ' + imagedir+figname+'n')

del ax
del fig

### ampere using B1 ###
for cp in ['x', 'y', 'z']:
    fig, axes = plt.subplots(figsize=(12,4), nrows=1, ncols=2, dpi=300)

    axes[0].plot(distance, data['curl_B1bats_'+cp], '+', color='DarkRed', label='$curl(B1)$')
    axes[0].plot(distance, unitmu0*data['Jbats_'+cp], 'x', color='Orange', label='$\\mu_0 J$')
    axes[0].set_yscale('symlog', linthreshy=1e-2)
    axes[0].legend()
    axes[0].set_title('dimensionfull values')
    axes[0].set_xlabel('distance from center [$R_E$]')
    axes[0].set_ylabel('values in  $\\frac{nT}{R_E}$')

    axes[1].plot(distance, (data['curl_B1bats_'+cp] - unitmu0*data['Jbats_'+cp])/(unitmu0*normJ),'x',
                 label='$\\frac{curl(B1) - \mu_0 J}{\mu_0 norm(J)}$', color='LightBlue')
    axes[1].plot(distance, (data['curl_B1bats_'+cp] - unitmu0*data['Jbats_'+cp])/(unitmu0*data['Jbats_'+cp]),'+',
                 label='$\\frac{curl(B1) - \mu_0 J}{\mu_0 J}$', color='Orange')
    axes[1].set_yscale('symlog', linthreshy=1e-4)
    axes[1].legend()
    axes[1].set_title('fractional errors')
    axes[1].set_xlabel('distance from center [$R_E$]')

    fig.suptitle(cp+' component', fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    figname = 'curlB1_and_J_'+cp+'_component.png'
    fig.savefig(imagedir+figname)
    if debug: print('saved png ' + imagedir+figname)
    log.write('saved png ' + imagedir+figname+'n')

    del axes
    del fig

amp_err = np.sqrt( (data['curl_B1bats_x'] - unitmu0*data['Jbats_x'])**2
                  +(data['curl_B1bats_y'] - unitmu0*data['Jbats_y'])**2
                  +(data['curl_B1bats_z'] - unitmu0*data['Jbats_z'])**2 )

amp_percent_err = 100.*amp_err/(unitmu0*normJ)

fig, ax = plt.subplots(figsize=(12,4), dpi=300)

ax.plot(distance, amp_percent_err, '.')
ax.set_yscale('symlog', linthreshy=1e-2)
ax.vlines(rCurrents, 0, 100, linewidths=1.)
ax.set_title('percent error in amperes law using only B1')
ax.set_xlabel('distance from center [$R_E$]')
ax.set_ylabel('percent error (%)')

figname = 'curlB1_and_J_percent_error.png'
fig.savefig(imagedir+figname)
if debug: print('saved png ' + imagedir+figname)
log.write('saved png ' + imagedir+figname+'n')

del ax
del fig



now = datetime.now()
log.write('script ended ' + now.strftime("%Y-%m-%d T%H:%M:%S") + '\n')
log.close()
