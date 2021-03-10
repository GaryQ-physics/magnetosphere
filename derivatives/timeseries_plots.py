import os
import sys
import numpy as np
import pandas as pd
import datetime

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util

run = 'DIPTSUR2'
n_times = 6


epsilons = [1./16., 1./8., 1./4., 1./2., 1., 2., 4., 8.]
df_variables = ['b1x', 'b1y', 'b1z', 'bx', 'by', 'bz', 'curl_b1_x', 'curl_b1_y',
       'curl_b1_z', 'curl_b_x', 'curl_b_y', 'curl_b_z', 'curl_j_x',
       'curl_j_y', 'curl_j_z', 'div_b', 'div_b1', 'div_j', 'div_jR',
       'gridspacing', 'jRx', 'jRy', 'jRz', 'jx', 'jy', 'jz', 'x', 'y',
       'z', 'jR_error']

times = list(util.get_available_slices(run)[1])
if n_times is not None:
    times = times[:n_times]

columns = []
for var in df_variables:
    columns.append('%s_sum_tot'%(var))
    columns.append('%s_sum_abs_tot'%(var))
    columns.append('%s_sum_sqr_tot'%(var))
    columns.append('%s_num_tot'%(var))
    for epsilon in epsilons:
        columns.append('%s_sum_epsilon_%f'%(var, epsilon))
        columns.append('%s_sum_abs_epsilon_%f'%(var, epsilon))
        columns.append('%s_sum_sqr_epsilon_%f'%(var, epsilon))
        columns.append('%s_num_epsilon_%f'%(var, epsilon))
timeseries = pd.DataFrame(columns=columns, index=range(len(times)))

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
    #        meta[key] = int(value)

    df = pd.read_pickle(fname_df)
    df['jR_error'] = np.sqrt( (df['jRx']-df['jx'])**2 \
                            + (df['jRy']-df['jy'])**2 \
                            + (df['jRz']-df['jz'])**2   )#!!!

    for var in df_variables:
        tmp = df[var]
        tmp = tmp[~np.isnan(tmp)]
        timeseries['%s_sum_tot'%(var)][i] = np.sum(tmp)
        timeseries['%s_sum_abs_tot'%(var)][i] = np.sum(np.abs(tmp))
        timeseries['%s_sum_sqr_tot'%(var)][i] = np.sum(tmp**2)
        timeseries['%s_num_tot'%(var)][i] = tmp.size
        for epsilon in epsilons:
            tmp = df[var][df['gridspacing'] == epsilon]
            tmp = tmp[~np.isnan(tmp)]
            timeseries['%s_sum_epsilon_%f'%(var, epsilon)][i] = np.sum(tmp)
            timeseries['%s_sum_abs_epsilon_%f'%(var, epsilon)][i] = np.sum(np.abs(tmp))
            timeseries['%s_sum_sqr_epsilon_%f'%(var,epsilon)][i] = np.sum(tmp**2)
            timeseries['%s_num_epsilon_%f'%(var, epsilon)][i] = tmp.size

timeseries.to_pickle('timeseries.pkl')
