"""
Generate cut plane for each file and for each variable associate with BATSRUS
run in conf['run_path'] specified in config.py. Files are downloaded as-needed.

Use cutplane_plot_demo2_animate.py to create .mp4 from images.

First set of cut planes at low resolution is used to compute min/max of
each variable. Second set is at high resolution with fixed axes limits
with limits determined using low resolution min/max values.

For long runs, execute from command line - Spyder often freezes after 50+ files
are processed.

When running in parallel wih Spyder/IPython, sometimes there is a broken pipe
error message. Program always runs after second attempt, however. In addition,
the disabiling of output buffering will not work (see line that starts with
sys.stdout = ...). 

TODO: If number of cores > number of vars, parallelize by files instead of
variables?

TODO: Generalize to use OpenGGCM and LFM. Main modification needed is the
generalization of CDF file names.
"""

import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from niceticks import niceticks
import cutplane_plot as cp
from util import filelist, filename2time, dlfile, time2datetime, filemeta
from events import events

event_list = events()

debug         = False  # Show extra logging information
para          = True   # Process in parallel using all CPUs
showplot      = False  # Show plot on screen. Set to false for long runs.
                       # Does not always work in Spyder/IPython, especially when
                       # para=True. Starting a new console sometimes fixes.
                       # TODO: Read PNG and display using PIL instead.
plot_type     = 2

# Testing options
first_only    = False  # Do only low-res first processing
second_only   = False  # Execute only high-res second processing
test_serial   = True  # Process two files in serial
test_parallel = False  # Process two files in parallel

vars = ['bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p','e']

opts = {
        'regen': True,       # Regenerate image even if found
        'plane': 'xz',       # Cut plane to plot
        'xlims': [-30, 15],  # Horizontal axis limits in R_E
        'ylims': [-20, 20],  # Vertical axis limits in R_E
        'zticks': None,      # z-axis ticks for variable
        'nf': None,          # Number of files to process. None => all files
        'dpi1': 96,          # Make multiple of 16 (for mimwrite animation)
        'dpi2': 96*3,        # Make multiple of 16 (for mimwrite animation)
        'delta1': 0.5,       # Low-res cut plane resolution in R_E
        'delta2': 0.1,       # High-res cut plane resolution in R_E
        'showplot': showplot # Show the plot on screen
        }

if first_only and second_only:
    raise ValueError("Both first_only and second_only are True. Only one may be True.")

event_list[0,0:5] = [2003,11,20,7,1]

if test_serial:
    para = False
    vars = ['p']
    opts["nf"] = 10
    opts["showplot"] = True

if test_parallel:
    para = True
    vars = ['p', 'jy']
    opts["nf"] = 3
    opts["showplot"] = True
    #opts["nf"] = 50
    #opts["showplot"] = False
    

def process_var(var, opts):

    k = 1
    times = []
    minmax = {'min': np.nan*np.empty(len(files)),
              'max': np.nan*np.empty(len(files)),
              'probe': {
                          'ux': np.nan*np.empty(len(files)),
                          'bz': np.nan*np.empty(len(files))
                    }
              }

    if plot_type == 2:
        dto = time2datetime(filename2time(files[0]))
        dtf = time2datetime(filename2time(files[-1]))
        dts = []
        data = {'ux': [], 'bz': []}

    for filename in files:

        if k > opts['nf']:
            break

        filename_png = conf["run_path_derived"] + "cutplanes/" + var + "/" \
                        + '{0:s}-{1:s}-type_{2:d}_delta_{3:.3f}.png' \
                        .format(filename, var, plot_type, opts['delta'])

        pkl = conf["run_path_derived"] + "cutplanes/minmax/" + var + '.pkl'

        zticks = opts['zticks']
        if os.path.exists(pkl):
            with open(pkl, 'rb') as handle:
                minmax = pickle.load(handle)

            if zticks is None:
                var_min = np.nanmin(minmax['min'])
                var_max = np.nanmax(minmax['max'])
                if not np.isnan(var_min) and not np.isnan(var_max):
                    zticks = niceticks(var_min, var_max, 8, debug=debug)

        if not opts['regen'] and os.path.exists(filename_png):
            # Image already exists
            print('Image exists for file {0:d}/{1:d} for {2:s} at resoluton {3:.2f}. Will not regenerate.'\
                  .format(k, opts['nf'], var, opts['delta']))
            k = k + 1
            continue

        times.append(filename2time(filename))

        dlfile(times[-1], debug=debug)

        if plot_type == 1:
            info = cp.plot(times[-1], var, opts['plane'],
                             dx=opts['delta'], dy=opts['delta'],
                             xlims=opts['xlims'], ylims=opts['ylims'],
                             dpi=opts['dpi'], showplot=opts['showplot'],
                             zticks=zticks,
                             png=True, pngfile=filename_png, debug=debug)
        elif plot_type == 2:

            def set_ylims(p, axes=None):
                if not np.all(np.isnan(minmax['probe'][p])):
                    var_min = np.nanmin(minmax['probe'][p])
                    var_max = np.nanmax(minmax['probe'][p])
                    if var_min != var_max:
                        yticks = niceticks(var_min, var_max, 8, debug=debug)
                        axes.set_yticks(yticks)
                        axes.set_ylim([var_min, var_max])

            # Not sure why next two code lines needed. Is it due to pandas import
            # in another hapiclient function? Issue this addresses described
            # at https://www.gitmemory.com/issue/facebook/prophet/999/500035319
            from pandas.plotting import register_matplotlib_converters
            register_matplotlib_converters()
            
            from matplotlib import pyplot as plt
            from hapiclient.plot.datetick import datetick
            from probe import probe
        
            fig = plt.figure(figsize=(7, 9), dpi=opts['dpi'])
            axes = fig.add_axes([0.07, 0.55, 0.9, 0.4])
            info = cp.plot(times[-1], var, opts['plane'],
                             axes=axes, 
                             dx=opts['delta'], dy=opts['delta'],
                             xlims=opts['xlims'], ylims=opts['ylims'],
                             logz=True,
                             showplot=opts['showplot'],
                             zticks=zticks,
                             debug=debug)
            
            time = filename2time(filename)
            dts.append(time2datetime(time))
            
            #import pdb;pdb.set_trace()
            # Sample upstream solar wind
            # TODO: Generalize code to allow plotting arbitrary number
            # of variables.
            d = probe(time, (31.5, 0, 0), ['ux', 'bz'])

            data['ux'].append(d['ux'])
            data['bz'].append(d['bz'])

            meta = filemeta(filename)
        
            axes = fig.add_axes([0.07, 0.1, 0.9, 0.18])
            axes.plot(dts, data['bz'])
            axes.set_xlim([dto, dtf])
            axes.grid()
            name = meta["parameters"]['bz']['plot_name']
            units = meta["parameters"]['bz']['vis_unit']
            axes.legend([name + ' [' + units + ']'])
            datetick('x', axes=axes, debug=False)
            set_ylims('bz', axes=axes)
            if np.all(times[-1][0:5] == event_list[0,0:5]):
                ylims = axes.get_ylim()
                axes.plot([dts[-1], dts[-1]], ylims, 'k-')
                
            axes = fig.add_axes([0.07, 0.28, 0.9, 0.18])
            axes.plot(dts, data['ux'])
            axes.grid()
            name = meta["parameters"]['ux']['plot_name']
            units = meta["parameters"]['ux']['vis_unit']
            axes.legend([name + ' [' + units + ']'])
            axes.set_xlim([dto, dtf])
            datetick('x', axes=axes)
            axes.axes.xaxis.set_ticklabels([])
            set_ylims('ux', axes=axes)

                            
            plt.savefig(filename_png, dpi=opts['dpi']) 
            
            if opts['showplot']:
                # These often won't cause plot to show.
                plt.show()
                fig.canvas.draw()
                fig.show()
            else:
                plt.close()

        else:
            raise ValueError("plot_type = {0:d} is not valid".format(plot_type))

        # This information will be the same for all vars but will
        # be stored in each variable's .pkl file.
        minmax['probe']['ux'][k-1] = d['ux']
        minmax['probe']['bz'][k-1] = d['bz']

        minmax['min'][k-1] = info['min']
        minmax['max'][k-1] = info['max']

        if debug:
            print("Writing " + pkl)
        with open(pkl, 'wb') as handle:
            pickle.dump(minmax, handle)
        if debug:
            print("Wrote " + pkl)

        print('Finished file {0:d}/{1:d} for {2:s} at resoluton {3:.2f}'\
               .format(k, opts['nf'], var, opts['delta']))
        
        k = k + 1

    return minmax

def process_all(opts):
    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is not None and num_cores > len(vars):
            num_cores = len(vars)
        print('Parallel processing {0:d} variable(s) using {1:d} cores'\
              .format(len(opts.keys()), num_cores))
        results = Parallel(n_jobs=num_cores)(\
                    delayed(process_var)(var, opts[var]) for var in vars)
    else:
        print('Serial processing {0:d} variable(s).'.format(len(opts.keys())))
        i = 0
        results = []
        for var in vars:
            results.append(process_var(var, opts[var]))
            i = i + 1
    return results

# Disable output buffering (sometimes causes problems with parallel runs)
# https://stackoverflow.com/questions/107705/disable-output-buffering
try:
    # Throws error in IPython, so put in try/catch.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
except:
    pass

# Create directory for variable min/max pkl files if needed
filename_path = conf["run_path_derived"] + "cutplanes/minmax/"
if not os.path.exists(filename_path):
    os.makedirs(filename_path)

# Create output directory if needed
filename_path = conf["run_path_derived"] + "cutplanes/"
if not os.path.exists(filename_path):
    os.makedirs(filename_path)
print('Image files will be written to subdirectory of\n' + filename_path)
# Create directory for each variable
for var in vars:
    filename_path = conf["run_path_derived"] + "cutplanes/" + var
    if not os.path.exists(filename_path):
        os.makedirs(filename_path)

files = filelist()

print('{0:d} files found.'.format(len(files)))

if opts['nf'] is None:
    # Process all files
    opts['nf'] = len(files)
else:
    print('Processing only first {0:d} files.'.format(opts['nf']))


# Create options dict for each variable.
Opts = {}
for var in vars:
    # .copy to ensure thread safety.
    # Currently not needed b/c opts is not modified.
    Opts[var] = opts.copy()

if not second_only:
    print('First processing.')
    for var in vars:
        Opts[var]['delta'] = opts['delta1']
        Opts[var]['dpi'] = opts['dpi1']

    process_all(Opts)

if not first_only:
    print('Second processing.')        
    for var in vars:
        Opts[var]['delta'] = opts['delta2']
        Opts[var]['dpi'] = opts['dpi2']

    process_all(Opts)
