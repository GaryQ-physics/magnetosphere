import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import util
import cxtransform as cx
import read_mag_grid_files as rmg
import magnetometers as mag

from hapiclient.plot.datetick import datetick
import matplotlib.pyplot as plt
import datetime


def main(run, LAT, LON , NStep=20, debug=False, plot_kwargs=None):
    direct = f'/tmp/{run}_{LAT}_{LON}/'
    if not os.path.exists(direct): os.mkdir(direct)

    times_fname = direct + 'times_for_dB.txt'
    if not os.path.exists(times_fname):
        #times1 = util.get_available_station_times(run, station)
        #filenames_all, times2 = util.get_available_slices(run)

        #t1 = times1[:,5] + 60*times1[:,4] + 60*60*times1[:,3] + 60*60*24*times1[:,2]  # ignores fraction of seconds !!!!
        #t2 = times2[:,5] + 60*times2[:,4] + 60*60*times2[:,3] + 60*60*24*times2[:,2]

        #t_common, ind1, ind2 = np.intersect1d(t1, t2, return_indices=True)
        #cdf_filenames = filenames_all[ind2]

        #assert(np.all(times2[ind2, 0:6] == times1[ind1, 0:6])) #!!!ignoring miliseconds!!!
        #time_common = times2[ind2, 0:6]

        cdf_filenames, time_common = util.get_available_slices(run)

        time_common = time_common[::NStep, :]
        cdf_filenames = cdf_filenames[::NStep]
        assert(time_common.shape[0] == cdf_filenames.size)

        with open(times_fname, 'w') as f:
            f.write('new_run:' + run + '\n')
            f.write('-- times\n')
            f.write('year, day, month, hr, min, sec\n')
            np.savetxt(f, time_common, fmt='%d')
    else:
        time_common = np.loadtxt(times_fname, skiprows=3)


    magfile_fname = direct + 'dB_from_magfile.txt'
    if not os.path.exists(magfile_fname):
        dB_magfile = np.nan*np.empty((time_common.shape[0], 12))
        for i in range(time_common.shape[0]):
            time = time_common[i,:]
            data, headers, csyst = rmg.getdata(run,time)
            kIndex, __spot_on = rmg.find_index(data, headers, LAT, LON)
            dB_magfile[i, :] = data[kIndex, 5:]
            assert(headers[5:] == ( 'dBnMhd', 'dBeMhd', 'dBdMhd',
                                    'dBnFac', 'dBeFac', 'dBdFac',
                                    'dBnHal', 'dBeHal', 'dBdHal',
                                    'dBnPed', 'dBePed', 'dBdPed') )

        with open(magfile_fname, 'w') as f:
            f.write('new_run:' + run + '\n')
            f.write('-- from mag file\n')
            f.write(', '.join(headers[5:])+'\n')
            np.savetxt(f, dB_magfile)

        print(np.all(dB_magfile[:,6:] == 0.))
    else:
        dB_magfile = np.loadtxt(magfile_fname, skiprows=3)
        headers, csyst = rmg.getdata(run,time_common[0,:], meta_only=True)

    calc_fname = direct + 'dB_calculated.txt'
    if not os.path.exists(calc_fname):
        dB_calc = np.nan*np.empty((time_common.shape[0], 12))
        for i in range(time_common.shape[0]):
            time = time_common[i,:]

            if (LAT,LON)==(18.907, 72.815) and csyst=='GEO': #colaba
                namefromstitch = conf[run+'_derived']+'timeseries/slices/' \
                        + 'B_biotsavart_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
                        + '_obs_point=colaba_rcut=%f.npy'%(util.get_rCurrents(run))
                dB_calc[i, 0:3] = (8**3/(8-2)**3) * np.load(namefromstitch) #!!!!!!!!! to account for the fact that calculated integral excluded the cells at block edges

            #dB_calc[i, 0:3] = 
            #dB_calc[i, 3:6] = 
            #dB_calc[i, 6:9] = 
            #dB_calc[i, 9:12] = 

        with open(calc_fname, 'w') as f:
            f.write('new_run: ' + run + '\n')
            f.write('-- calculated\n')
            f.write('dBnMHD, dBeMHD, dBdMHD, dBnFAC, dBeFAC, dBdFAC, ...\n')
            np.savetxt(f, dB_calc)
    else:
        dB_calc = np.loadtxt(calc_fname, skiprows=3)

    ########## plotting #############
    if plot_kwargs is None:
        plot_kwargs = {'component': 'norm', 'contributor':'mhd'}

    dtimes = []
    for i in range(time_common.shape[0]):
        dtimes.append(datetime.datetime(*time_common[i,0:6]))

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
    plt.title(f"deltaB at {LAT},{LON}  {plot_kwargs['contributor']} {plot_kwargs['component']}")
    plt.legend()
    plt.savefig(f'{direct}dB_plot.png')
    #plt.show()
    plt.clf()

if __name__ == '__main__':
    main('DIPTSUR2', 18.907, 72.815)
