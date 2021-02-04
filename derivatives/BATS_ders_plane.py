import os
import numpy as np
import pandas as pd

from config import conf
from probe import GetRunData
import util
from units_and_constants import phys

from derivatives import GetDivergence, GetCurl, GetFrobeniusNormDel, GetOperatorNormDel


####################
run = 'DIPTSUR2'
#pntlist = 'native_random_sampled2'
plane = 'xz'
offset = 1./16.
epsilon = 1./8.
para = True
debug = False
####################

pntlist = '%s_plane_y=%f'%(plane,offset)

if run == 'DIPTSUR2':
    #time = (2019,9,2,6,30,0,0)
    time = (2019,9,2,4,10,0,0)
    rCurrents = 1.8
    rBody = 1.5

direct = conf[run+'_derived'] + 'derivatives/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6)
direct = direct + pntlist + '/'

points_fname = conf['storage'] + 'pointlists/' + pntlist + '__points.txt'
points = np.loadtxt(points_fname)
if debug: print(points)

def GetDerivativesArray(epsil):
    results_fname = direct + 'partial_derivatives_epsilon=%f.bin'%(epsil)

    results = np.fromfile(results_fname).reshape((3, points.shape[0], 3, 3))
    if debug: print('loaded "results" array from ' + results_fname + '\n')

    pointsloaded = np.fromfile(direct + 'derivatives_points.bin').reshape((points.shape[0], 3))
    assert(np.all(points == pointsloaded))

    return results


###### process data #############
dels = [GetDerivativesArray(epsilon)]

del_j_batsrus = dels[0][0,:,:,:]
del_b_batsrus = dels[0][1,:,:,:]
del_b1_batsrus = dels[0][2,:,:,:]

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

distance = np.sqrt(data['point_x']**2 + data['point_y']**2 + data['point_z']**2)
if debug:
    print('distance:')
    print(distance)

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


ng = 42
#### 2d color plot #####
data = ON_div_B1bats.values.reshape((ng,ng))
var = 'B1'
###########

print(type(data))
data = np.nan_to_num(data)

x_grid = points[:, 0].reshape((ng,ng))
z_grid = points[:, 2].reshape((ng,ng))

import matplotlib.pyplot as plt

fig = plt.figure()
axes = fig.gca()

#axes.set_title(title, fontsize=10, family='monospace')

Nbt = 4 # Approximate number of colors between ticks for linear scale
import matplotlib
cmap = matplotlib.pyplot.get_cmap('viridis', Nbt*(1000-1))

axes.set(title='color plot %s_y=%f epsilon=%f'%(var,offset,epsilon))
axes.set(xlabel=plane[0])
axes.set(ylabel=plane[1])
axes.axis('square')
axes.set_xlim(-8., 8.)
axes.set_ylim(-8., 8.)

datamax = np.max(data)
datamin = np.min(data)
print(datamax,datamin)

import matplotlib.colors as colors
norm = colors.SymLogNorm(linthresh=0.01, vmin=datamin, vmax=datamax)
pcm = axes.pcolormesh(x_grid, z_grid, data, norm=norm, cmap=cmap, vmin=datamin, vmax=datamax)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(axes)
cax = divider.append_axes("right", size="5%", pad=0.05)

cb = axes.figure.colorbar(pcm, cax=cax, label='Normalized Divergence %s'%(var))

imagedir = conf['base'] + 'images/' + run + '/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6)
fig.savefig(imagedir+'color_plot_%s_y=%f_epsilon=%f.png'%(var,offset,epsilon))

#plt.show()
