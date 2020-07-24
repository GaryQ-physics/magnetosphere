#SCARR5
import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import biot_savart_kameleon_interpolated_grid as bsk
import util
import cxtransform as cx
import read_mag_grid_files as rmg


samp = 50*np.arange(17)
#YKClat = 62.480
#YKClon = 245.518
MLAT = 0.
MLON = 0.

def commonTimes(debug=False):
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
    #times2 = np.array(util.timelist(listtxt='ls-1.txt'))

    t1 = times1[:,5] + 60*times1[:,4] + 60*60*times1[:,3] + 60*60*24*times1[:,2]  # ignores fraction of seconds !!!!
    t2 = times2[:,5] + 60*times2[:,4] + 60*60*times2[:,3] + 60*60*24*times2[:,2]
    t_common, ind1, ind2 = np.intersect1d(t1, t2, return_indices=True)

    dB_SWMF = np.nan*np.empty((samp.size, 3))
    for i in range(samp.size):
        dB_SWMF[i, :] = np.array(rmg.analyzedata( 
                                            str(magfilenames[ind1[samp[i]]]),
                                            MLAT, MLON, debug=debug) )[1:4]


    if debug:
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

    assert(np.all(times1[ind1, 0:6] == times2[ind2, 0:6])) #!!!ignoring miliseconds!!!
    #assert(tuple(headers[13:16]) == ('dBn', 'dBe', 'dBd'))
    filenames_all = np.array(filenames_all, dtype=str)
    return [times1[ind1, :], filenames_all[ind2], dB_SWMF]  # dBn dBe dBd  as  0,1,2


def dB_kam_tofile(time_common, filenames, debug=False, tag=None, xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.125):
    assert(time_common.shape[0] == filenames.shape[0])
    if debug:
        print(time_common)
        print(filenames)
        print(dB_SWMF)

    if tag == None:
        tup = xlims+ylims+zlims+(d,)
        tag = '_{0:07.2f}_{1:07.2f}_{2:07.2f}_{3:07.2f}_{4:07.2f}_{5:07.2f}_{6:.5f}_'.format(*tup)

    datafname = conf['run_path_derived'] +'dB_kam_tofile' + tag + '.txt'
    f = open(datafname,'a') # append only mode
    f.write('\n  new_run: xlims=' + str(xlims) + ', ylims=' + str(ylims) + ', zlims=' + str(zlims) + ', d=' + str(d) + '\n')
    f.write('-- $' + str(d) + 'R_E$ interpolated grid; (X,Y,Z) = (' + str(xlims) + ', ' + str(ylims) + ', ' + str(zlims) + ')\n')
    f.write('year, day, month, hr, min, sec, dB_kam_0, dB_kam_1, dB_kam_2\n')

    dB_kam = np.nan*np.empty((samp.size, 3))
    for i in range(samp.size):
        print('i,samp[i] = ' + str((i,samp[i])) )
        time = time_common[samp[i],:]
        filename = str(filenames[samp[i]]) # numpy string --> normal string

        util.dlfile(filename, debug=True)

        #mpos = cx.GEOtoMAG([1., YKClat, YKClon] , time, 'sph', 'sph')
        #mlat = mpos[1]
        #mlon = mpos[2]

        print(time)
        print(filename)

        dB_kam[i,:] = bsk.run(time, MLAT, MLON, filename=conf['run_path'] + filename, para=True,
            xlims=xlims, ylims=ylims, zlims=zlims, d=d,
            print_output=True)

        f.write(str(time)[1:-1] + '   ' + str(dB_kam[i,:])[1:-1] + '\n')

    f.close()
    return dB_kam

def dB_kam_fromfile(filename, fullname=False):
    if not fullname:
        filename = conf['run_path_derived'] + filename
    arr = np.loadtxt(filename, skiprows=4)
    times = arr[:, 0:6]
    dB_kam = arr[:, 6:9]
    f = open(filename, 'r')
    f.readline()
    f.readline()
    label = f.readline()
    return [times, dB_kam, label]


def plot_from_file(datafnames, fullname=False, component='norm', compsyst='MAG'):
    from hapiclient.plot.datetick import datetick

    import matplotlib.pyplot as plt
    import datetime

    if type(datafnames) == str:
        datafnames = [datafnames]

    time_common, filenames, dB_SWMF = commonTimes()
    #dB_SWMF = dB_SWMF_allcommon[samp, :]

    for datafname in datafnames:
        times, dB_kam, label = dB_kam_fromfile(datafname, fullname=fullname)

        tup = tuple([float(datafname.split('_')[i+3]) for i in range(7)])

        assert(times.shape[0] == 17)
        assert(dB_kam.shape[0] == 17)

        times = np.array(times, dtype=int)

        station_pos = cx.MAGtoGSM(np.array([1., MLAT, MLON]), times, 'sph', 'car')
        if compsyst == 'MAG':
            Pole = cx.MAGtoGSM(np.array([0., 0., 1.]), times, 'car', 'car')
        elif compsyst == 'GEO':
            Pole = cx.GEOtoGSM(np.array([0., 0., 1.]), times, 'car', 'car')
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
            dB_kam_comp = np.sqrt(np.einsum('ij,ij->i', dB_kam, dB_kam))
        else:
            if component == 'east':
                u = U1
            elif component == 'north':
                u = U2
            elif component == 'down':
                u = -U3
            dB_kam_comp = np.einsum('ij,ij->i', u, dB_kam)

        plt.plot(dtimes, dB_kam_comp, label=label, alpha=0.7)


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





    datetick('x')
    # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
    #plt.gcf().autofmt_xdate()

    plt.xlabel('time')
    plt.xlabel('dB (in nanoTesla)')
    plt.title('dB ' + component + ' using ' + compsyst + ' system')
    plt.legend()
    plt.show()




