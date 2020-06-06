"""
Generate cut plane for each file for each variable.

First set of cut planes  at low resolution is used to compute min/max of
each variable. Second set is at high resolution. Use cut_plane_plot_animate.py
to create movie from images.

When running in parallel wih IPython, sometimes thre is a broken pipe error
message. Program always runs after second attempt, however.
"""

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from niceticks import niceticks
import cut_plane_plot as cp
from urlretrieve import urlretrieve

###############################################################################
para = True   # Process in parallel using all CPUs
debug = False
vars = ['bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p','e']
###############################################################################
opts = {
        'regen': True,      # Regenerate image even if found
        'plane': 'xz',      # Cut plane to plot
        'xlims': [-30, 15], # Horizontal axis limits in R_E
        'ylims': [-20, 20], # Vertical axis limits in R_E
        'zticks': None,     # z-axis ticks for variable
        'dpi': 300,
        'nf': None,         # Number of files to process. None => all files
        'delta1': 0.5,      # Low-res cut plane resolution in R_E
        'delta2': 0.1,      # High-res cut plane resolution in R_E
        }

def filename2ts(filename):
    # Extract time stamp from file name
    tstr = filename[11:] 
    y, m, d = int(tstr[0:4]), int(tstr[4:6]), int(tstr[6:8])
    h, M, s = int(tstr[9:11]), int(tstr[11:13]), int(tstr[13:15])
    f = int(tstr[16:19])
    return [y, m, d, h, M, s, f]

def process_var(var, opts):

    k = 1
    times = []
    maxes = []
    mins = []
    for file in files:
        if k > opts['nf']:
            break
        filename = file.rstrip()
        if not filename[0:9] == '3d__var_3':
            continue
        k = k + 1

        filename_png = filename_path + filename + '-' + var \
                        + '-' + '{0:.3f}'.format(opts['delta']) + '.png'

        if not opts['regen'] and os.path.exists(filename_png):
            # Image already exists
            continue

        times.append(filename2ts(filename))

        fname_full = conf['run_path'] + filename
        if not os.path.exists(fname_full):
            fileurl = conf['run_url'] + filename
            if debug:
                print('Downloading ' + fileurl)
                print('to')
            fname_tmp = fname_full + ".tmp"
            if debug:
                print(fname_tmp)
            urlretrieve(fileurl, fname_tmp)
            if debug:
                print('Downloaded ' + fileurl)
            os.rename(fname_tmp, fname_full)
            if debug:
                print('Renamed *.cdf.tmp to *.cdf')

        # TODO: Catch download error
        if os.path.exists(fname_full):
            info = cp.plot(times[-1], opts['plane'], None, var, 
                             dx=opts['delta'], dy=opts['delta'],
                             xlims=opts['xlims'], ylims=opts['ylims'],
                             dpi=opts['dpi'],
                             zticks=opts['zticks'],
                             png=True, pngfile=filename_png, debug=debug)

        mins.append(info['min'])
        maxes.append(info['max'])
        print('Finished file {0:03d}/{1:d} for {2:s} at resoluton {3:.2f}'\
               .format(k, opts['nf'], var, opts['delta']))
    return mins, maxes

def process_all(opts):
    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        print('Parallel processing {0:d} variables using {1:d} cores'\
              .format(len(opts.keys()), num_cores))
        results = Parallel(n_jobs=num_cores)(\
                    delayed(process_var)(var, opts[var]) for var in opts.keys())
    else:
        print('Serial processing {0:d} variables.'.format(len(opts.keys())))
        i = 0
        results = []
        for var in opts.keys():
            results.append(process_var(var, opts[var]))
            i = i + 1
    return results

# Create output directory if needed
filename_path = conf["run_path_derived"] + "cutplanes/"
if not os.path.exists(filename_path):
    os.makedirs(filename_path)
print('Image files will be written to\n' + filename_path)

# Get list of run files
urlretrieve(conf['run_url'] + 'ls-1.txt', conf['run_path'] + 'ls-1.txt')

# Read list of run files
with open(conf['run_path'] + 'ls-1.txt','r') as f:
    files = f.readlines()

if opts['nf'] is None:
    # Process all files
    opts['nf'] = len(files)

opts['delta'] = opts['delta1']
Opts = {} 
for var in vars:
    Opts[var] = opts

print('First processing.')
results = process_all(Opts)

results = np.array(results)
ticks = []
i = 0
for var in vars:
    var_min = np.min(results[i,:][0])
    var_max = np.max(results[i,:][1])
    ticks = niceticks(var_min, var_max, 8, debug=False)
    Opts[var]['zticks'] = ticks
    Opts[var]['delta'] = opts['delta2']

print('Second processing.')
results = process_all(Opts)
