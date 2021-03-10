import numpy as np
import pandas as pd
from native_grid import main

import util

run = 'DIPTSUR2'
n_times = 500
para = True

times = list(util.get_available_slices(run)[1])
if n_times is not None:
    times = times[:n_times]

if para:
    from joblib import Parallel, delayed
    import multiprocessing
    num_cores = multiprocessing.cpu_count()
    if num_cores is not None and num_cores > len(times):
        num_cores = len(times)
    print('Parallel processing {0:d} time slices using {1:d} cores'\
          .format(len(times), num_cores))
    timeseries = Parallel(n_jobs=num_cores)(\
                    delayed(main)(run, time) for time in times)
else:
    timeseries = []
    for time in times:
        timeseries.append(main(run, time))


#print(timeseries)
print(len(timeseries[0].keys()))

timeDF = pd.DataFrame(columns=timeseries[0].keys(), index=range(len(times)))

for key in timeDF.keys():
    for i in range(len(times)):
        timeDF[key][i] = timeseries[i][key]

timeDF.to_pickle('timeDF.pkl')
