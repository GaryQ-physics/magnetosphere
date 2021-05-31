import os
import sys
import numpy as np
from hapiclient.plot.datetick import datetick
import matplotlib.pyplot as plt
import datetime

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
import read_mag_grid_files as rmg
import read_computed as rcp


def plot_magfile_calc(run, surface_location , NStep=20, debug=False, plot_kwargs=None):
    direct = conf[run+'_derived'] + f'rmf_{surface_location}/'
    if not os.path.exists(direct): os.mkdir(direct)

    cdf_filenames, times = util.get_available_slices(run)
    times = times[::NStep, :]
    cdf_filenames = cdf_filenames[::NStep]

    dB_magfile = np.empty((times.shape[0],12))
    dB_calc = np.empty((times.shape[0],12))

    for i in range(times.shape[0]):
        time = times[i,:]
        dB_magfile[i, :] = rmg.get_mag_grid_values(run, time, surface_location)

        dB_calc[i, 0:3] = rcp.get_mhd_values(run, time, surface_location , "B_biotsavart")
        dB_calc[i, 3:6] = rcp.get_mhd_values(run, time, surface_location , "B_fac", norcut=True)
        #dB_calc[i, 6:9] = 
        #dB_calc[i, 9:12] = 


    ########## plotting #############
    if plot_kwargs is None:
        plot_kwargs = {'component': 'norm', 'contributor':'mhd'}


    dtimes = []
    for i in range(times.shape[0]):
        dtimes.append(datetime.datetime(*times[i,0:6]))

    if plot_kwargs['contributor'] == 'mhd':
        offset = 0
    elif plot_kwargs['contributor'] == 'fac':
        offset = 3

    for i in range(2):
        if i==0:
            label = 'magfile'
            dB = dB_magfile
        elif i==1:
            label = 'calculated'
            dB = dB_calc

        if plot_kwargs['component'] == 'norm':
            dB_comp = np.sqrt(np.einsum('ij,ij->i', dB[:, offset:offset+3], dB[:, offset:offset+3]))
        else:
            if plot_kwargs['component'] == 'north':
                dB_comp = dB[:, offset+0]
            elif plot_kwargs['component'] == 'east':
                dB_comp = dB[:, offset+1]
            elif plot_kwargs['component'] == 'down':
                dB_comp = dB[:, offset+2]

        plt.plot(dtimes, dB_comp, label=label, alpha=0.7)

    #datetick('x')
    # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
    #plt.gcf().autofmt_xdate()
    plt.xlabel('time')
    plt.ylabel('nano Tesla')
    plt.title(f"deltaB at {surface_location}  {plot_kwargs['contributor']} {plot_kwargs['component']}")
    plt.legend()
    plt.savefig(f"{direct}dB_plot_{plot_kwargs['contributor']}_{plot_kwargs['component']}.png")
    #plt.show()
    plt.clf()

if __name__ == '__main__':
    plot_magfile_calc('DIPTSUR2', "colaba", plot_kwargs={'component': 'norm', 'contributor':'mhd'})
    plot_magfile_calc('DIPTSUR2', "colaba", plot_kwargs={'component': 'norm', 'contributor':'fac'})
