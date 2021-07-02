import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd

#from hapiclient.plot.datetick import datetick
import magnetopost.util as mputil

import read_mag_grid_files as rmg

def plot_magfile_calc(rundir, surface_location, component='norm'):
    run = mputil.prep_run(rundir)

    magnetotimes = run['magnetosphere_files'].keys()
    ionotimes = run['ionosphere_files'].keys()

    dB_magfile = np.empty((len(magnetotimes),12))
    cachefile = f'/tmp/dB_magfile-{surface_location}.npy'

    if os.path.exists(cachefile):
        dB_magfile = np.load(cachefile)
    else:
        indexes = np.empty(len(magnetotimes), dtype=np.int32)
        for i,time in enumerate(magnetotimes):
            filename = run['magnetosphere_files'][time]
            filename = filename.replace('3d__var_2_e','mag_grid_e')[:-8]+'.out'
            headers, csyst = rmg.getmetadata(None, None, filename=filename)
            data = rmg.getdata(None, None, filename=filename)
            assert(headers[5:] == ( 'dBnMhd', 'dBeMhd', 'dBdMhd',
                                    'dBnFac', 'dBeFac', 'dBdFac',
                                    'dBnHal', 'dBeHal', 'dBdHal',
                                    'dBnPed', 'dBePed', 'dBdPed') )
            _r, LAT, LON = mputil.GetMagnetometerCoordinates(surface_location, time, csyst, 'sph')
            if LON < 0:
                LON = LON + 360.
            if csyst != 'GEO':
                print(csyst)
            if abs(LAT-18.907)>1e-8 or LON != 72.815:
                print('warning')
                print(LAT,LON)

            kIndex, __spot_on = rmg.find_index(data, headers, LAT, LON)
            dB_magfile[i, :] = data[kIndex, 5:]
            indexes[i] = kIndex

        np.save(cachefile, dB_magfile)
        np.save(f'/tmp/indexes-{surface_location}.npy', indexes)

    d_magnetotimes = []
    i = -1
    for time in magnetotimes:
        i += 1
        d_magnetotimes.append(datetime.datetime(*time))

    d_ionotimes = []
    i = -1
    for time in ionotimes:
        i += 1
        d_ionotimes.append(datetime.datetime(*time))

    def get_comp(dB, comp):
        if comp == 'norm':
            dB_comp = np.sqrt(np.einsum('ij,ij->i', dB, dB))
        else:
            print('WARNING not implemented properly yet')
            if comp == 'north':
                dB_comp = dB[:, 0]
            elif comp == 'east':
                dB_comp = dB[:, 1]
            elif comp == 'down':
                dB_comp = dB[:, 2]
        return dB_comp
    

    fig, axes = plt.subplots(figsize=(10, 5), nrows=4, ncols=1, dpi=200, sharex=True)
    offset = -3
    axnum = -1
    for contributor, dtimes in [('msph', d_magnetotimes), ('fac', d_magnetotimes), ('hall', d_ionotimes), ('pedersen', d_ionotimes)]:
        offset += 3
        axnum += 1

        #missfact = 1.
        #if axnum == 3:
        #    missfact = 1./(4.*np.pi)

        axes[axnum].plot(d_magnetotimes, get_comp(
                                    dB_magfile[:, offset:offset+3],
                                                component), label='magfile', alpha=0.7)
        axes[axnum].plot(dtimes, get_comp(
                np.load(f'{rundir}/derived/timeseries/bs_{contributor}-{surface_location}.npy'),
                                         component), label='calculated', alpha=0.7)

        axes[axnum].set_ylabel('nano Tesla')
        axes[axnum].set_title(f"contibution of {contributor}")

    axes[-1].set_xlabel('time')

    fig.tight_layout(rect=[0,0.03,1,0.95])
    #datetick('x')
    # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
    #plt.gcf().autofmt_xdate()
    fig.suptitle(f"deltaB at {surface_location}. Plotting {component}")
    fig.legend()
    fig.savefig(f"{rundir}/derived/dB_plot_{surface_location}_{component}.png")
    fig.savefig(f"{rundir}/derived/dB_plot_{surface_location}_{component}.pdf")
    #plt.show()
    plt.close(fig)


def onlymsph_plot_magfile_calc(rundir, surface_location, component='norm'):
    run = mputil.prep_run(rundir)

    magnetotimes = run['magnetosphere_files'].keys()

    dB_magfile = np.empty((len(magnetotimes),12))
    cachefile = f'/tmp/dB_magfile-{surface_location}.npy'

    if os.path.exists(cachefile):
        dB_magfile = np.load(cachefile)
    else:
        indexes = np.empty(len(magnetotimes), dtype=np.int32)
        for i,time in enumerate(magnetotimes):
            filename = run['magnetosphere_files'][time]
            filename = filename.replace('3d__var_2_e','mag_grid_e')[:-8]+'.out'
            headers, csyst = rmg.getmetadata(None, None, filename=filename)
            data = rmg.getdata(None, None, filename=filename)
            assert(headers[5:] == ( 'dBnMhd', 'dBeMhd', 'dBdMhd',
                                    'dBnFac', 'dBeFac', 'dBdFac',
                                    'dBnHal', 'dBeHal', 'dBdHal',
                                    'dBnPed', 'dBePed', 'dBdPed') )
            _r, LAT, LON = mputil.GetMagnetometerCoordinates(surface_location, time, csyst, 'sph')
            if LON < 0:
                LON = LON + 360.
            if csyst != 'GEO':
                print(csyst)
            if abs(LAT-18.907)>1e-8 or LON != 72.815:
                print('warning')
                print(LAT,LON)

            kIndex, __spot_on = rmg.find_index(data, headers, LAT, LON)
            dB_magfile[i, :] = data[kIndex, 5:]
            indexes[i] = kIndex

        np.save(cachefile, dB_magfile)
        np.save(f'/tmp/indexes-{surface_location}.npy', indexes)

    d_magnetotimes = []
    i = -1
    for time in magnetotimes:
        i += 1
        d_magnetotimes.append(datetime.datetime(*time))

    def get_comp(dB, comp):
        if comp == 'norm':
            dB_comp = np.sqrt(np.einsum('ij,ij->i', dB, dB))
        else:
            print('WARNING not implemented properly yet')
            if comp == 'north':
                dB_comp = dB[:, 0]
            elif comp == 'east':
                dB_comp = dB[:, 1]
            elif comp == 'down':
                dB_comp = dB[:, 2]
        return dB_comp
    

    fig, axes = plt.subplots(figsize=(10, 5), nrows=2, ncols=1, dpi=200, sharex=True)
    contributor = 'msph'
    dtimes = d_magnetotimes
    offset = 0

    calc = get_comp(
            np.load(f'{rundir}/derived/timeseries/bs_{contributor}-{surface_location}.npy'), component)
    mag = get_comp(dB_magfile[:, offset:offset+3], component)

    axes[0].plot(d_magnetotimes, mag, label="SWMF magfile's dBMhd contribution",
            color='LightBlue', marker='x', alpha=1, linestyle='None', markersize=3)
    axes[0].plot(dtimes, calc, label=r'Calculated |$\int_{\mathcal{M}}$biotsavart|',
            color='Orange', marker='+', alpha=1, linestyle='None', markersize=3)
    axes[0].set_ylabel('nT')
    axes[0].legend()
    axes[0].set_title(f"Magnetosphere Contribution to $|B(x_0)|$ where x_0 is Colaba")

    axes[1].plot(dtimes, calc-mag, label='difference')
    axes[1].set_ylabel('nT')
    axes[1].set_xlabel('time')
    axes[1].legend()
    axes[1].set_title(f"Difference of above plots (calc-magfile)")

    fig.tight_layout(rect=[0,0.03,1,0.95])
    #datetick('x')
    # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
    #plt.gcf().autofmt_xdate()
    fig.savefig(f"{rundir}/derived/dBmhd_plot_{surface_location}_{component}.png")
    fig.savefig(f"{rundir}/derived/dBmhd_plot_{surface_location}_{component}.pdf")
    #plt.show()
    plt.close(fig)


if __name__ == '__main__':
    pth = '/home/gary/media_sunspot/'
    onlymsph_plot_magfile_calc(pth+'DIPTSUR2/', "colaba", component='norm')
    onlymsph_plot_magfile_calc(pth+'DIPTSUR2/', "gridpnt1", component='norm')
    onlymsph_plot_magfile_calc(pth+'DIPTSUR2/', "gridpnt2", component='norm')
    #plot_magfile_calc(pth+'DIPTSUR2/', "colaba")
