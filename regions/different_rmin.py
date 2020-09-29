import os
import sys
import numpy as np
import pickle

import regions
import magnetometers as mg

#https://stackoverflow.com/questions/16779497/how-to-set-memory-limit-for-thread-or-process-in-python
import resource
soft, hard = 72*2**30, 72*2**30
resource.setrlimit(resource.RLIMIT_AS,(soft, hard))


para = True
serial = False
run = 'IMP10_RUN_SAMPLE'
time = (2019,9,2,7,0,0)
location = mg.GetMagnetometerLocation('colaba', (2019,1,1,1,0,0), 'MAG', 'sph')

rs = 0.03125*np.arange(32,36)

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

pkl = 'different_rmin.pkl'

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





#toret[0,2,:]   dBs[:,0,2,:] is Nx3

