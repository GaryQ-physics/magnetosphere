import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import biot_savart_kameleon as bsk
import util
import cxtransform as cx
import read_ccmc_datafiles as r_ccmc

samp = 50*np.arange(31)
YKClat = 62.480
YKClon = 245.518

def commonTimes(debug=False):
    data, headers = r_ccmc.getdata(2006)

    listnames = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/SWPC_SWMF_052811_2_GM_cdf_list'
    a = np.loadtxt(listnames, dtype=str, skiprows=1) #(_,5)

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


def dB_kam_tofile(time_common, filenames, debug=False, tag=None, xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.125):
    assert(time_common.shape[0] == filenames.shape[0])
    if debug:
        print(time_common)
        print(filenames)
        print(dB_SWMF)

    if tag == None:
        tup = xlims+ylims+zlims+(d,)
        tag = '_{0:07.2f}_{1:07.2f}_{2:07.2f}_{3:07.2f}_{4:07.2f}_{5:07.2f}_{6:.5f}_'.format(*tup)

    datafname = conf['data_path'] + 'dB_kam_tofile' + tag + '.txt'
    f = open(datafname,'a') # append only mode
    f.write('\n  new_run: xlims=' + str(xlims) + ', ylims=' + str(ylims) + ', zlims=' + str(zlims) + ', d=' + str(d) + '\n')
    f.write('-- $' + str(d) + 'R_E$ interpolated grid; (X,Y,Z) = (' + str(xlims) + ', ' + str(ylims) + ', ' + str(zlims) + ')')
    f.write('year, day, month, hr, min, sec, dB_kam_0, dB_kam_1, dB_kam_2\n')

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
            xlims=xlims, ylims=ylims, zlims=zlims, d=d,
            print_output=True)

        f.write(str(time)[1:-1] + '   ' + str(dB_kam[i,:])[1:-1] + '\n')

    f.close()
    return dB_kam

def dB_kam_fromfile(filename, fullname=False):
    if not fullname:
        filename = conf['data_path'] + filename
    arr = np.loadtxt(filename, skiprows=4)
    times = arr[:, 0:6]
    dB_kam = arr[:, 6:9]
    f = open(filename, 'r')
    f.readline()
    f.readline()
    label = f.readline()
    return [times, dB_kam, label]


def plot_from_file(datafnames, fullname=False):
    import matplotlib.pyplot as plt
    #import matplotlib.dates
    import datetime

    if type(datafnames) == str:
        datafnames = [datafnames]

    time_common, filenames, dB_SWMF_allcommon = commonTimes()

    dB_SWMF = dB_SWMF_allcommon[samp, :]

    for datafname in datafnames:
        times, dB_kam, label = dB_kam_fromfile(datafname, fullname=fullname)

        tup = tuple([float(datafname.split('_')[i+3]) for i in range(7)])

        assert(times.shape[0] == 31)
        assert(dB_kam.shape[0] == 31)

        times = np.array(times, dtype=int)
        dtimes = []
        for i in range(times.shape[0]):
            dtimes.append(datetime.datetime(times[i,0],times[i,1],times[i,2],times[i,3],times[i,4],times[i,5]))

        dB_kam_norm = np.sqrt(np.einsum('ij,ij->i', dB_kam, dB_kam))
        plt.plot(dtimes, dB_kam_norm, label=label)

    times = np.array(time_common[samp,:], dtype=int)
    dtimes = []
    for i in range(times.shape[0]):
        dtimes.append(datetime.datetime(times[i,0],times[i,1],times[i,2],times[i,3],times[i,4],times[i,5]))

    dB_SWMF_norm = np.sqrt(np.einsum('ij,ij->i', dB_SWMF, dB_SWMF))
    plt.plot(dtimes, dB_SWMF_norm, label='-- SWMF')


    # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
    plt.gcf().autofmt_xdate()

    plt.xlabel('year = 2006')
    plt.legend()
    plt.show()

'''
import time as tm

to = tm.time()

time_common, filenames, dB_SWMF_allcommon = commonTimes()
ret = dB_kam_tofile(time_common, filenames, tag=None, xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.25)

tf = tm.time()

print(ret)
print('time = ' + str((tf-to)/3600.) + ' hours') # time = 3.4527500342 hours
'''

plot_from_file(['dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.50000_.txt',
            'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.25000_.txt',
            'dB_kam_tofile_-056.00_0008.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt'])
#       or alternatively
#plot_from_file('/home/gary/magnetosphere/dB_kam_tofile_copy.txt', fullname=True)
