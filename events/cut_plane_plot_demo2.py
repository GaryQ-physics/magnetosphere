"""
Generate cut plane for each file for each variable.

First set of cut planes  at low resolution is used to compute min/max of
each variable. Second set is at high resolution. Use cut_plane_plot_animate.py
to create movie from images.

When running in parallel wih IPython, sometimes thre is a broken pipe error
message. Program always runs after second attempt, however. In addition,
the disabiling of output buffering will not work (see line that starts with
sys.stdout =).

TODO: If number of cores > number of vars, parallelize by files instead of
variables.
"""

import os
import sys
import numpy as np
import pickle

#from pympler import tracker
#tr = tracker.SummaryTracker()

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from niceticks import niceticks
import cut_plane_plot as cp
from urlretrieve import urlretrieve

# https://stackoverflow.com/questions/107705/disable-output-buffering
# Throws error in IPython.
try:
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
except:
    pass

debug = False
para = True   # Process in parallel using all CPUs

vars = ['bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p','e']
if False:
    num_cores = None
    vars = ['p']

opts = {
        'regen': False,     # Regenerate image even if found
        'plane': 'xz',      # Cut plane to plot
        'xlims': [-30, 15], # Horizontal axis limits in R_E
        'ylims': [-20, 20], # Vertical axis limits in R_E
        'zticks': None,     # z-axis ticks for variable
        'nf': 50,           # Number of files to process. None => all files
        'dpi1': 112,        # Make multiple of 16 (for mimwrite)
        'dpi2': 112*3,        # Make multiple of 16 (for mimwrite)
        'delta1': 0.5,      # Low-res cut plane resolution in R_E
        'delta2': 0.1,     # High-res cut plane resolution in R_E
        'showplot': True   # Show the plot on screen
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

        filename_png = conf["run_path_derived"] + "cutplanes/" + var + "/" \
                        + filename + '-' + var \
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
                             dpi=opts['dpi'], showplot=opts['showplot'],
                             zticks=opts['zticks'],
                             png=True, pngfile=filename_png, debug=debug)

        mins.append(info['min'])
        maxes.append(info['max'])
        print('Finished file {0:d}/{1:d} for {2:s} at resoluton {3:.2f}'\
               .format(k, opts['nf'], var, opts['delta']))
        
        k = k + 1

    return mins, maxes

def process_all(opts):
    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is not None and num_cores > len(vars):
            num_cores = len(vars)
        print('Parallel processing {0:d} variables using {1:d} cores'\
              .format(len(opts.keys()), num_cores))
        results = Parallel(n_jobs=num_cores)(\
                    delayed(process_var)(var, opts[var]) for var in vars)
    else:
        print('Serial processing {0:d} variables.'.format(len(opts.keys())))
        i = 0
        results = []
        for var in vars:
            results.append(process_var(var, opts[var]))
            i = i + 1
    return results

# Create output directory if needed
filename_path = conf["run_path_derived"] + "cutplanes/"
if not os.path.exists(filename_path):
    os.makedirs(filename_path)
print('Image files will be written to\n' + filename_path)
for var in vars:
    filename_path = conf["run_path_derived"] + "cutplanes/" + var
    if not os.path.exists(filename_path):
        os.makedirs(filename_path)

# Get list of run files
urlretrieve(conf['run_url'] + 'ls-1.txt', conf['run_path'] + 'ls-1.txt')

# Read list of run files
with open(conf['run_path'] + 'ls-1.txt','r') as f:
    files = f.readlines()

if opts['nf'] is None:
    # Process all files
    opts['nf'] = len(files)

opts['delta'] = opts['delta1']
opts['dpi'] = opts['dpi1']
Opts = {} 
for var in vars:
    Opts[var] = opts.copy()

pkl = 'cut_plane_plot_demo2.pkl'
if '__file__' in globals():
    pkl = os.path.dirname(os.path.abspath(__file__)) + '/' + pkl

if True:
    print('First processing.')
    results = process_all(Opts)
    
    # Save results file in case error occurs after first processing.
    with open(pkl, 'wb') as handle:
        pickle.dump(results, handle)

if True:
    print('Second processing.')
    with open(pkl, 'rb') as handle:
        results = pickle.load(handle)
    
    from matplotlib import pyplot as plt
    results = np.array(results)
    ticks = []
    i = 0
    Opts.clear()
    for var in vars:
        var_min = np.min(results[i,:][0])
        var_max = np.max(results[i,:][1])
        if False:
            plt.close()
            plt.title(var)
            plt.hist(results[i,:][0], color='b', edgecolor='b')
            plt.hist(results[i,:][1], color='r', edgecolor='r')
            plt.legend({'min','max'})
            plt.show()
        i = i + 1
        ticks = niceticks(var_min, var_max, 8, debug=False)
        Opts[var] = opts.copy()
        Opts[var]['zticks'] = ticks.copy()
        Opts[var]['delta'] = opts['delta2']
        Opts[var]['dpi'] = opts['dpi2']

    process_all(Opts)
