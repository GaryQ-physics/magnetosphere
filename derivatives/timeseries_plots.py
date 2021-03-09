import os
import sys
import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util

run = 'DIPTSUR2'
n_times = None

epsilons = [1./16., 1./8., 1./4., 1./2., 1., 2., 4., 8.]

times = list(util.get_available_slices(run)[1])
if n_times is not None:
    times = times[:n_times]

keys = []
keys.append('integral')
keys.append('absintegral')
keys.append('num_total')
for epsilon in epsilons:
    keys.append('sum_epsilon_%f'%(epsilon))
    keys.append('abssum_epsilon_%f'%(epsilon))
    keys.append('num_epsilon_%f'%(epsilon))

timeseries = pd.DataFrame(columns=keys, index=range(len(times)))

for i in range(len(times)):
    if i%10 == 0: print('i=%d'%(i))
    time = times[i]
    direct = conf[run+'_derived'] + 'derivatives/native_grid/'
    fname_df = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_df.pkl'%util.tpad(time, length=6)
    fname_meta = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_meta.txt'%util.tpad(time, length=6)

    #with open(fname_meta,'r') as handle:
    #    meta = {}
    #    for line in handle.readlines():
    #        key, value = line.split(' ')
    #        meta[key]=int(value)

    df = pd.read_pickle(fname_df)

    tmp = df['div_b1']*df['gridspacing']
    tmp = tmp[tmp!=np.nan]
    timeseries['integral'][i] = np.sum(tmp)
    timeseries['absintegral'][i] = np.sum(np.abs(tmp))
    timeseries['num_total'][i] = tmp.size

    for epsilon in epsilons:
        tmp = df['div_b1'][df['gridspacing'] == epsilon]
        tmp = tmp[tmp!=np.nan]
        timeseries['sum_epsilon_%f'%(epsilon)][i] = np.sum(tmp)
        timeseries['abssum_epsilon_%f'%(epsilon)][i] = np.sum(np.abs(tmp))
        timeseries['num_epsilon_%f'%(epsilon)][i] = tmp.size

timeseries.to_pickle('timeseries.pkl')



dtimes = []
for time in times:
    dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))

