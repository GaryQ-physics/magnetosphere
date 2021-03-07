from native_grid_nan_ends import main

import util

run = 'DIPTSUR2'
n_times = 2
para = False

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
    Parallel(n_jobs=num_cores)(\
            delayed(main)(run, time) for time in times)
else:
    for time in times:
        main(run, time)

