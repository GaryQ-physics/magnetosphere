import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import util
import regions
import cxtransform as cx
import magnetometers as mg

run = 'IMP10_RUN_SAMPLE'
pm = 16.125
reg =  {'xlims': (-pm, pm),
        'ylims': (-pm, pm),
        'zlims': (-pm, pm),
        'd': 0.25
        }
para = True
serial = False
rs = 0.03125*np.arange(32,64) #96

pkl = run + '_different_rmin_pm_%f_d_%f.pkl'%(pm, reg['d'])


if run == 'DIPTSUR2':
    time = (2019,9,2,6,30,0,0)
    SWMF = np.array([-1058. , 2. ,  -22.])
if run == 'IMP10_RUN_SAMPLE':
    time = (2019,9,2,7,0,0,0)
    SWMF = np.array([-1021., 2.24, 121.69])

location = mg.GetMagnetometerLocation('colaba', (2019,1,1,1,0,0), 'MAG', 'sph')

pkl = conf[run+'_derived'] \
      + 'regions/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6) \
      + pkl

if para:
    from joblib import Parallel, delayed
    import multiprocessing
    num_cores = multiprocessing.cpu_count()
    if num_cores is not None and num_cores > len(rs):
        num_cores = len(rs)
    print('Parallel processing {0:d} file(s) using {1:d} cores'\
          .format(len(rs), num_cores))
    dBs = Parallel(n_jobs=num_cores)(\
            delayed(regions.signedintegrate)(run, time, location,\
            regions=(reg,), fwrite=False, rmin=r) for r in list(rs))
elif serial:
    print('Serial processing {0:d} file(s)'\
          .format(len(rs)))
    dBs = []
    for r in list(rs):
        dBs.append(regions.signedintegrate(run, time, location, regions=(reg,), fwrite=False, rmin=r))

if para or serial:
    dBs = np.array(dBs)
    print("Writing " + pkl)
    with open(pkl, 'wb') as handle:
        pickle.dump(dBs, handle)

if (not para) and (not serial):
    if os.path.exists(pkl):
        with open(pkl, 'rb') as handle:
            print("Reading " + pkl)
            dBs = pickle.load(handle)
    else:
        raise ValueError('Cannot find needed pickle file')

print(dBs[:,0,2,:])
print(dBs.shape)

calc = dBs[:,0,2,:]
swmf = np.repeat([SWMF], dBs.shape[0], axis=0)

normdiff =  calc - swmf
normdiff = np.sqrt(np.einsum('ij,ij->i', normdiff, normdiff))

diffnorm = np.sqrt(np.einsum('ij,ij->i', calc, calc)) - np.sqrt(np.einsum('ij,ij->i', swmf, swmf))

title = run + '_%04d-%02d-%02dT%02d:%02d:%02d.%03d' % tuple(time)
title = title + '\nat mlat=%.3f, mlon=%.3f, MLT=%.3f hrs'%(location[1], location[2], cx.MAGtoMLT(location[2], time))


if os.path.exists(conf[run + '_cdf'] + '../../PARAM.in'):
    f = open(conf[run + '_cdf'] + '../../PARAM.in', 'r')
    lines = f.readlines()
    Rcurrents = None
    Rbody = None
    for line in lines:
        if 'Rcurrents' in line:
            Rcurrents = float(line.split(' ')[0])
        if 'Rbody' in line:
            Rbody = float(line.split(' ')[0])
else:
    Rcurrents = None
    Rbody = None

try:
    import matplotlib.pyplot as plt

    plt.plot(rs, normdiff)
    plt.hlines(1., 0., 3., colors='k')
    if Rcurrents is not None:
        plt.vlines(Rcurrents, 0., 450., label='Rcurrents', colors='r')
    if Rbody is not None:
        plt.vlines(Rbody, 0., 450., label='Rbody', colors='g')
    plt.legend()
    plt.title(title)
    plt.xlabel('rmin $R_E$')
    plt.ylabel('$|dB_{calculated} - dB_{swmf}|$ [nT]')
    print('saving ' + pkl + '-normdiff.png')
    plt.savefig(pkl + '-normdiff.png')
    plt.clf()

    plt.plot(rs, diffnorm)
    plt.hlines(1., 0., 3., colors='k')
    if Rcurrents is not None:
        plt.vlines(Rcurrents, 0., 450., label='Rcurrents', colors='r')
    if Rbody is not None:
        plt.vlines(Rbody, 0., 450., label='Rbody', colors='g')
    plt.legend()
    plt.title(title)
    plt.xlabel('rmin $R_E$')
    plt.ylabel('$|dB_{calculated}| - |dB_{swmf}|$ [nT]')
    print('saving ' + pkl + '-diffnorm.png')
    plt.savefig(pkl + '-diffnorm.png')
    plt.clf()
except:
    print('Failed to plot. Can run again with para = False and serial = False to skip computation')


