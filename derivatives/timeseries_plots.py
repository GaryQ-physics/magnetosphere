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
epsilons = [1./16., 1./8., 1./4., 1./2.]

times = list(util.get_available_slices(run)[1])
if n_times is not None:
    times = times[:n_times]

averages = np.empty((len(epsilons), len(times))); averages[:,:]=np.nan
integrals = np.empty(len(times)); integrals[:]=np.nan
dtimes = []
for i in range(len(times)):
    time = times[i]
    dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
    direct = conf[run+'_derived'] + 'derivatives/native_grid/'
    fname_df = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_df.pkl'%util.tpad(time, length=6)
    fname_meta = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_meta.txt'%util.tpad(time, length=6)

    with open(fname_meta,'r') as handle:
        meta = {}
        for line in handle.readlines():
            key, value = line.split(' ')
            meta[key]=int(value)

    df = pd.read_pickle(fname_df)

    tmparr = np.array(df['div_b1']*df['gridspacing'])
    tmparr = np.nan_to_num(tmparr)
    tmparr = np.abs(tmparr)
    integrals[i] = np.sum(tmparr)

    for j in range(len(epsilons)):
        tmp = df['div_b1'][df['gridspacing'] == epsilons[j]]
        tmp = tmp[tmp!=np.nan]
        summ = np.sum(tmp)
        assert(len(tmp.shape)==1)
        nn = tmp.size
        averages[j,i] = summ/nn


plt.plot(dtimes, integrals,'.',)
plt.savefig('integrals.png')
plf.clf()

plt.plot(dtimes, averages[0,:],'.',)
plt.savefig('averages0.png')
plf.clf()

plt.plot(dtimes, averages[1,:],'.',)
plt.savefig('averages1.png')
plf.clf()

plt.plot(dtimes, averages[2,:],'.',)
plt.savefig('averages2.png')
plf.clf()

plt.plot(dtimes, averages[3,:],'.',)
plt.savefig('averages3.png')
plf.clf()
