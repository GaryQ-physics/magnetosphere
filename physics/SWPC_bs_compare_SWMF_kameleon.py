import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import biot_savart_kameleon as bsk
import util
import cxtransform as cx
import read_ccmc_datafiles as r_ccmc


def commonTimes(debug=False):
    data, headers = r_ccmc.getdata(2006)

    listnames = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/SWPC_SWMF_052811_2_GM_cdf_list'
    a = np.loadtxt(listnames, dtype=str, skiprows=1) #(N,5)

    times1 = np.array(data[:, 0:6], dtype=int)
    times2 = np.column_stack([np.array(list(np.char.split(a[:,2],sep='/')), 
        dtype=int), np.array(list(np.char.split(a[:,4],sep=':')), dtype=int)])

    t1 = times1[:,5] + 60*times1[:,4] + 60*60*times1[:,3] + 60*60*24*times1[:,2]
    t2 = times2[:,5] + 60*times2[:,4] + 60*60*times2[:,3] + 60*60*24*times2[:,2]
    t_common, ind1, ind2 = np.intersect1d(t1, t2, return_indices=True)

    if debug:
        print(a.shape)
        print(t1.shape)
        print(t1)
        print(t2.shape)
        print(t2)
        print(t_common.shape)
        print(t_common)
        print(ind1.shape)
        print(ind1)
        print(ind2.shape)
        print(ind2)
        print(times1[ind1, :])
        print(times2[ind2, :])

    assert(np.all(times1[ind1, :] == times2[ind2, :]))
    assert(tuple(headers[13:16]) == ('dBn', 'dBe', 'dBd'))
    return [times1[ind1, :], a[ind2,0], data[ind1, 13:16]]


def dB_kam_tofile(time_common, filenames, debug=False, tag=''):
    assert(time_common.shape[0] == filenames.shape[0])
    if debug:
        print(time_common)
        print(filenames)
        print(dB_SWMF)

    datafname = conf['data_path'] + 'dB_kam_tofile' + tag + '.txt'
    f = open(datafname,'a') # append only mode
    f.write('\n  new_run: xlims=(-56., 8.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.125\n')
    f.write('time, dB_kam_0, dB_kam_1, dB_kam_2\n')

    samp = 50*np.arange(31)
    #samp = 500*np.arange(3)

    YKClat = 62.480
    YKClon = 245.518
    dB_kam = np.nan*np.empty((samp.size, 3))

    for i in range(samp.size):
        print('i,samp[i] = ' + str((i,samp[i])) )
        time = time_common[samp[i],:]
        filename = filenames[samp[i]]

        util.dlfile_SWPC(filename, debug=True)

        mpos = cx.GEOtoMAG([1., YKClat, YKClon] , time, 'sph', 'sph')
        mlat = mpos[1]
        mlon = mpos[2]

        print(time)
        print(filename)

        dB_kam[i,:] = bsk.run(time, mlat, mlon, filename='/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/' + filename, para=True,
            xlims=(-56., 8.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.125,
            print_output=True)

        f.write(str(time) + ' ')
        f.write(str(dB_kam[i,:]) + '\n')
        #f.write('%e %e %e'%(dB_kam[i,0], dB_kam[i,1], dB_kam[i,2]))

    f.close()
    return dB_kam

def dB_kam_fromfile(filename, fullname=False):
    if not fullname:
        filename = conf['data_path'] + filename
    arr = np.loadtxt(filename, skiprows=3)
    times = arr[:, 0:6]
    dB_kam = arr[:, 6:9]
    return [times, dB_kam]


def plot_from_file(filename):
    time, dB_kam = dB_kam_fromfile(filename)
    time_common, filenames, dB_SWMF_allcommon = commonTimes()

    samp = 50*np.arange(31)
    dB_SWMF = dB_SWMF_allcommon[samp, :]

    import matplotlib.pyplot as plt
    #import matplotlib.dates
    import datetime

    time = np.array(time, dtype=int)
    dtimes = []
    for i in range(time.shape[0]):
        dtimes.append(datetime.datetime(time[i,0],time[i,1],time[i,2],time[i,3],time[i,4],time[i,5]))

    #dtimes = matplotlib.dates.date2num(dtimes)

    dB_kam_norm = np.sqrt(np.einsum('ij,ij->i', dB_kam, dB_kam))
    dB_SWMF_norm = np.sqrt(np.einsum('ij,ij->i', dB_SWMF, dB_SWMF))
    plt.plot(dtimes, dB_kam_norm, label='dB_kam_norm')
    plt.plot(dtimes, dB_SWMF_norm, label='dB_SWMF_norm')


    # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
    plt.gcf().autofmt_xdate()

    plt.xlabel('year = 2006')
    plt.legend()
    plt.show()
    
'''
import time as tm

to = tm.time()

time_common, filenames, dB_SWMF_allcommon = commonTimes()
ret = dB_kam_tofile(time_common,filenames)

tf = tm.time()

print(ret)
print('time = ' + str((tf-to)/3600.) + ' hours') # time = 3.4527500342 hours
'''

plot_from_file('dB_kam_tofile.txt')
