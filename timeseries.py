'''man
 -s : only stitch, do not regenerate slices. (if slices not already on disk, this will error out)
 -p : compute in parallel
 -v : verbose, print details ect
 --nfiles : set equal to n_files (default: None, i.e. all files)

eg:
python timeseries.py -pv --nfiles=4
'''
from job_config import *

##### interperate command line inputs ##################################
import os
import sys

n_times = None
para = False
verbose = False
stitch_only = False
for arg in sys.argv:
    if arg[0:2] == '--':
        if '--nfiles' == arg[0:8]:
            n_times = int(arg[9:])
        else:
            raise ValueError ('unknown argument '+arg)
    elif arg[0] == '-':
        for char in arg[1:]:
            if char=='s':
                stitch_only = True
            if char=='p':
                para = True
            if char=='v':
                verbose = True
########################################################################

from config import conf
import util
from swmf_file_reader.read_swmf_files import read_all
from integrals import slice_B_biotsavart, stitch_B_biotsavart, slice_B_coulomb, stitch_B_coulomb
from stats_summary import slice_stats_summary, stitch_stats_summary
from probe_locations import slice_probe_b1, stitch_probe_b1
from cutplane_native import slice_xzplane, stitch_xzplane

# if rcut=None or not supplied in the job_config file, it defaults to rCurrents
try:
    if rcut is None:
        rcut = util.get_rCurrents(run)
except NameError:
    rcut = util.get_rCurrents(run)

# make needed directories if not exists.
if not os.path.exists(conf[run+'_derived']+'timeseries/slices/'):
    os.makedirs(conf[run+'_derived']+'timeseries/slices/')
    print('created directory '+conf[run+'_derived']+'timeseries/slices/')

# wrapper function, to run all functions desired by command line argument,
#    while only loading the .out once per time slice. 
def wrap(time):
    cache = read_all(util.time2CDFfilename(run,time)[:-8])

    if do_cutplane:
        slice_xzplane(run, time, rcut=rcut, cache=cache)

    if do_stats_summary:
        slice_stats_summary(run, time, rcut=rcut, cache=cache)

    if do_biotsavart_integral:
        for point in integral_points:
            slice_B_biotsavart(run, time, point, rcut=rcut, cache=cache)

    if do_coulomb_integral:
        for point in integral_points:
            slice_B_coulomb(run, time, point, rcut=rcut, cache=cache)

    if do_probing:
        for point in probe_points:
            slice_probe_b1(run, time, point, cache=cache)

# loop through each time slice, in parallel or in serial, and execute wrapper
times = list(util.get_available_slices(run)[1])
if n_times is not None:
    times = times[:1*n_times:1]

if not stitch_only: # typically always runs unless overriden by -s option
    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is None:
            assert(False)
        num_cores = min(num_cores, len(times), 20)
        print('Parallel processing {0:d} time slices using {1:d} cores'\
              .format(len(times), num_cores))
        Parallel(n_jobs=num_cores)(\
                  delayed(wrap)(time) for time in list(times))
    else:
        for time in list(times):
            wrap(time)

# stitch the written files together
if do_cutplane:
    stitch_xzplane(run, times, rcut=rcut)

if do_stats_summary:
    stitch_stats_summary(run, times, rcut=rcut)

if do_biotsavart_integral:
    for point in integral_points:
        stitch_B_biotsavart(run, times, point, rcut=rcut)

if do_coulomb_integral:
    for point in integral_points:
        stitch_B_coulomb(run, times, point, rcut=rcut)

if do_probing:
    for point in probe_points:
        stitch_probe_b1(run, times, point)
