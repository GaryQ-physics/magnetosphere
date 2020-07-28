import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import biot_savart_kameleon_interpolated_grid as bsk
import util
import cxtransform as cx
import read_ccmc_datafiles as r_ccmc
import read_mag_grid_files as rmg


samp_SWPC = 50*np.arange(31)
#YKClat = 62.480
#YKClon = 245.518

magnetometer_stations = {'YKClat':62.480, 'YKClon':245.518}

samp_SCARR5 = 50*np.arange(17)
#YKClat = 62.480
#YKClon = 245.518
#MLAT = 0.
#MLON = 0.

def commonTimes(run, station, debug=False, skip_SWMF=False):
    """
    run = 'SWPC',  station = 'YKC'

    run = 'SCARR5',  station = [MLAT, MLON]  where these are floats
    """

    if run == 'SWPC':
        assert(type(station) == str)
    elif run == 'SCARR5':
        assert(type(station) == list)

    if run == 'SWPC':
        samp = samp_SWPC
        data, headers = r_ccmc.getdata([2006, station])

        listnames = conf['SWPC_cdf_path'] + 'SWPC_SWMF_052811_2_GM_cdf_list'
        a = np.loadtxt(listnames, dtype=str, skiprows=1) #(_,5)

        times1 = np.array(data[:, 0:6], dtype=int)
        times2 = np.column_stack([np.array(list(np.char.split(a[:,2],sep='/')), 
            dtype=int), np.array(list(np.char.split(a[:,4],sep=':')), dtype=int)])
    elif run == 'SCARR5':
        samp = samp_SCARR5
        magfilenames = np.loadtxt(conf['run_path'] + 'ls-1_magfiles.txt', dtype=str)[:,0]
        times1 = []
        for file_name in list(magfilenames):
            times1.append(rmg.mag_grid_file2time(file_name))
        times1 = np.array(times1)

        times2 = []
        filenames_all = util.filelist(listtxt = 'ls-1.txt')
        for file_name in filenames_all:
            times2.append(util.filename2time(file_name))
        times2 = np.array(times2)

    t1 = times1[:,5] + 60*times1[:,4] + 60*60*times1[:,3] + 60*60*24*times1[:,2]  # ignores fraction of seconds !!!!
    t2 = times2[:,5] + 60*times2[:,4] + 60*60*times2[:,3] + 60*60*24*times2[:,2]
    t_common, ind1, ind2 = np.intersect1d(t1, t2, return_indices=True)

    assert(np.all(times1[ind1, 0:6] == times2[ind2, 0:6])) #!!!ignoring miliseconds!!!

    if run == 'SWPC':
        directory = conf['SWPC_derived'] + station + '/'

        filenames = a[ind2,0]

    elif run == 'SCARR5':
        directory = conf['run_path_derived'] + 'mpos:{0:.3f}:{1:.3f}/'.format(*tuple(station))

        filenames_all = np.array(filenames_all, dtype=str)
        filenames = filenames_all[ind2]

    if not os.path.exists(directory):
        os.makedirs(directory)
        print('Created directory ' + directory)

    if not skip_SWMF:
        if run == 'SWPC':
            assert(tuple(headers[13:16]) == ('dBn', 'dBe', 'dBd'))
            dB_SWMF_allcommon = data[ind1, 13:16] # dBn dBe dBd  as  0,1,2
            dB_SWMF = dB_SWMF_allcommon[samp, :]
        elif run == 'SCARR5':
            dB_SWMF = np.nan*np.empty((samp.size, 3))
            for i in range(samp.size):
                dB_SWMF[i, :] = np.array(rmg.analyzedata( 
                                                    str(magfilenames[ind1[samp[i]]]),
                                                    station[0], station[1], debug=debug) )[1:4]

        datafname = directory + 'dB_SWMF_tofile' + '.txt'
        f = open(datafname,'a') # append only mode
        f.write('\n  new_run:\n')
        f.write('-- SWMF\n')
        f.write('year, day, month, hr, min, sec, dB_SWMF_0, dB_SWMF_1, dB_SWMF_2\n')
        for i in range(samp.size):
            f.write(str(times1[ind1[samp[i]], :])[1:-1] + '   ' + str(dB_SWMF[i,:])[1:-1] + '\n')
        f.close()

    return [times1[ind1, :], filenames]


def dB_kam_tofile(run, station, time_common, filenames, debug=False, tag=None, xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.125):
    assert(time_common.shape[0] == filenames.shape[0])
    if debug:
        print(time_common)
        print(filenames)
        print(dB_SWMF)
    assert(run == 'SWPC' or run == 'SCARR5')

    if tag == None:
        tup = xlims+ylims+zlims+(d,)
        tag = '_{0:07.2f}_{1:07.2f}_{2:07.2f}_{3:07.2f}_{4:07.2f}_{5:07.2f}_{6:.5f}_'.format(*tup)

    if run == 'SWPC':
        samp = samp_SWPC
        datafname = conf['SWPC_derived'] + station + '/'
    elif run == 'SCARR5':
        samp = samp_SCARR5
        directory = conf['run_path_derived'] + 'mpos:{0:.3f}:{1:.3f}/'.format(*tuple(station))

    datafname = directory + 'dB_kam_tofile' + tag + '.txt'
    f = open(datafname,'a') # append only mode
    f.write('\n  new_run: xlims=' + str(xlims) + ', ylims=' + str(ylims) + ', zlims=' + str(zlims) + ', d=' + str(d) + '\n')
    f.write('-- $' + str(d) + 'R_E$ interpolated grid; (X,Y,Z) = (' + str(xlims) + ', ' + str(ylims) + ', ' + str(zlims) + ')\n')
    f.write('year, day, month, hr, min, sec, dB_kam_0, dB_kam_1, dB_kam_2\n')

    dB_kam = np.nan*np.empty((samp.size, 3))
    for i in range(samp.size):
        print('i,samp[i] = ' + str((i,samp[i])) )
        time = time_common[samp[i],:]
        filename = str(filenames[samp[i]]) # numpy string -> normal string

        if run == 'SWPC':
            util.dlfile_SWPC(filename, debug=True)
            mpos = cx.GEOtoMAG([1., magnetometer_stations[station + 'lat'], magnetometer_stations[station + 'lon']] , time, 'sph', 'sph')
            mlat = mpos[1]
            mlon = mpos[2]
            filename_full = conf['SWPC_cdf_path'] + filename
        elif run == 'SCARR5':
            util.dlfile(filename, debug=True)
            mlat = station[0]
            mlon = station[1]
            filename_full = conf['run_path'] + filename


        print(time)
        print(filename)

        dB_kam[i,:] = bsk.run(time, mlat, mlon, filename=filename_full, para=True,
            xlims=xlims, ylims=ylims, zlims=zlims, d=d,
            print_output=True)

        f.write(str(time)[1:-1] + '   ' + str(dB_kam[i,:])[1:-1] + '\n')

    f.close()
    return dB_kam


def dB_fromfile(run, station, filename, fullname=False):
    if not fullname:
        if run == 'SWPC':
            filename = conf['SWPC_derived'] + station + '/' + filename
        elif run == 'SCARR5':
            filename = conf['run_path_derived'] + 'mpos:{0:.3f}:{1:.3f}/'.format(*tuple(station)) + filename
    arr = np.loadtxt(filename, skiprows=4)
    times = arr[:, 0:6]
    dB = arr[:, 6:9]
    f = open(filename, 'r')
    f.readline()
    f.readline()
    label = f.readline()
    return [times, dB, label]


def plot_from_file(run, station, datafnames, fullname=False, component='norm', compsyst='MAG'):
    """
    compsyst = 'MAG'  (default)
    compsyst = 'GEO'
    """
    from hapiclient.plot.datetick import datetick

    import matplotlib.pyplot as plt
    import datetime

    if type(datafnames) == str:
        datafnames = [datafnames]

    #time_common, filenames, dB_SWMF = commonTimes(run)

    for datafname in datafnames:
        times, dB, label = dB_fromfile(run, station, datafname, fullname=fullname)

        #tup = tuple([float(datafname.split('_')[i+3]) for i in range(7)])


        times = np.array(times, dtype=int)
        if compsyst == 'MAG':
            Pole = cx.MAGtoGSM(np.array([0., 0., 1.]), times, 'car', 'car')
        elif compsyst == 'GEO':
            Pole = cx.GEOtoGSM(np.array([0., 0., 1.]), times, 'car', 'car')

        if run == 'SWPC':
            assert(times.shape[0] == 31)
            assert(dB.shape[0] == 31)

            station_pos = cx.GEOtoGSM(np.array([1., magnetometer_stations[station + 'lat'], magnetometer_stations[station + 'lon']]), times, 'sph', 'car')

        elif run == 'SCARR5':
            assert(times.shape[0] == 17)
            assert(dB.shape[0] == 17)

            station_pos = cx.MAGtoGSM(np.array([1., station[0], station[1]]), times, 'sph', 'car')

        U3 = station_pos
        U3norm = np.sqrt(U3[:,0]**2 + U3[:,1]**2 + U3[:,2]**2)#( not really needed since norm is 1 in these units where R_e=1, but just to be consistent in all units)
        divU3norm = 1./U3norm
        U3 = U3*divU3norm[:,np.newaxis]
        U1 = np.cross(Pole, U3)
        U1norm = np.sqrt(U1[:,0]**2 + U1[:,1]**2 + U1[:,2]**2)
        divU1norm = 1./U1norm
        #https://stackoverflow.com/questions/5795700/multiply-numpy-array-of-scalars-by-array-of-vectors
        U1 = U1*divU1norm[:,np.newaxis]
        U2 = np.cross(U3, U1)
        U = np.column_stack([U1, U2, U3])

        dtimes = []
        for i in range(times.shape[0]):
            dtimes.append(datetime.datetime(times[i,0],times[i,1],times[i,2],times[i,3],times[i,4],times[i,5]))

        if component == 'norm':
            dB_comp = np.sqrt(np.einsum('ij,ij->i', dB, dB))
        else:
            if 'SWMF' in datafname:
                if component == 'east':
                    dB_SWMF_comp = dB_SWMF[:, 1]
                elif component == 'north':
                    dB_SWMF_comp = dB_SWMF[:, 0]
                elif component == 'down':
                    dB_SWMF_comp = dB_SWMF[:, 2]

            else:
                if component == 'east':
                    u = U1
                elif component == 'north':
                    u = U2
                elif component == 'down':
                    u = -U3
                dB_comp = np.einsum('ij,ij->i', u, dB)

        plt.plot(dtimes, dB_comp, label=label, alpha=0.7)

    '''
    times = np.array(time_common[samp,:], dtype=int)
    dtimes = []
    for i in range(times.shape[0]):
        dtimes.append(datetime.datetime(times[i,0],times[i,1],times[i,2],times[i,3],times[i,4],times[i,5]))

    if component == 'norm':
        dB_SWMF_comp = np.sqrt(np.einsum('ij,ij->i', dB_SWMF, dB_SWMF))
    else: # dBn dBe dBd  as  0,1,2
        if component == 'east':
            dB_SWMF_comp = dB_SWMF[:, 1]
        elif component == 'north':
            dB_SWMF_comp = dB_SWMF[:, 0]
        elif component == 'down':
            dB_SWMF_comp = dB_SWMF[:, 2]

    plt.plot(dtimes, dB_SWMF_comp, label='-- SWMF', alpha=0.7)
    '''

    datetick('x')
    # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
    #plt.gcf().autofmt_xdate()

    plt.xlabel('time')
    plt.ylabel('dB (in nanoTesla)')
    plt.title('dB ' + component + ' using ' + compsyst + ' system')
    plt.legend()
    plt.show()


def compute(run, station, debug=False, tag=None, xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.125, skip_SWMF=False):
    import time as tm
    to = tm.time()

    time_common, filenames = commonTimes(run, station, skip_SWMF=skip_SWMF)
    ret = dB_kam_tofile(run, station, time_common, filenames, debug=debug, tag=tag, xlims=xlims, ylims=ylims, zlims=zlims, d=d)

    tf = tm.time()

    print('time = ' + str((tf-to)/3600.) + ' hours') # time = 3.4527500342 hours
    g = open('done_time.txt', 'a')
    g.write('\n time = ' + str((tf-to)/3600.) + ' hours\n')
    g.close()
    print(ret)
