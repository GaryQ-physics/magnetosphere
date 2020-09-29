import os
import sys
import numpy as np
import pickle

import regions
import cxtransform as cx
import magnetometers as mg

#https://stackoverflow.com/questions/16779497/how-to-set-memory-limit-for-thread-or-process-in-python
import resource
soft, hard = 72*2**30, 72*2**30
resource.setrlimit(resource.RLIMIT_AS,(soft, hard))

run = 'IMP10_RUN_SAMPLE'

pkl = run + 'different_rmin.pkl'
para = True
serial = False
rs = 0.03125*np.arange(32,96)

if run == 'DIPTSUR2':
    time = (2019,9,2,6,30,0,0)
if run == 'IMP10_RUN_SAMPLE':
    time = (2019,9,2,7,0,0,0)
location = mg.GetMagnetometerLocation('colaba', (2019,1,1,1,0,0), 'MAG', 'sph')

pm = 31.875
reg =  {'xlims': (-pm, pm),
        'ylims': (-pm, pm),
        'zlims': (-pm, pm),
        'd': 0.25
        }

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

if (not para) and (not serial) and os.path.exists(pkl):
    with open(pkl, 'rb') as handle:
        print("Reading " + pkl)
        dBs = pickle.load(handle)

print(dBs[:,0,2,:])
print(dBs.shape)

if run == 'DIPTSUR2':
    SWMF = np.array([-1058. , 2. ,  -22.])
if run == 'IMP10_RUN_SAMPLE':
    SWMF = np.array([-1021., 2.24, 121.69])

toplot = dBs[:,0,2,:] - np.repeat([SWMF], dBs.shape[0], axis=0)
toplot = np.sqrt(np.einsum('ij,ij->i', toplot, toplot))
title = pkl[:-4] + '\n%04d-%02d-%02dT%02d:%02d:%02d.%03d' % tuple(time)
title = title + '\nat mlat=%.3f, mlon=%.3f, MLT=%.3f hrs'%(location[1], location[2], cx.MAGtoMLT(location[2], time))

import matplotlib.pyplot as plt
plt.plot(rs, toplot)
plt.title(title)
plt.xlable('rmin $R_E$')
plt.ylable('|dB_calculated - dB_swmf| $nT$')
print('saving ' +pkl + '.png')
plt.savefig(pkl + '.png')
plt.clf()



#toret[0,2,:]   dBs[:,0,2,:] is Nx3

