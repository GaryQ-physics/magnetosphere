import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import SWPC_bs_compare_SWMF_kameleon as swpc
#import cxtransform as cx
import time as tm


plot = True

def plot_explore():
    from hapiclient.plot.datetick import datetick
    import matplotlib.pyplot as plt
    import datetime

    times, dB_kam1, label1 = swpc.dB_kam_fromfile('dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.50000_.txt', fullname=False)
    assert(times.shape[0] == 31)
    times, dB_kam2, label2 = swpc.dB_kam_fromfile('dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.25000_.txt', fullname=False)
    assert(times.shape[0] == 31)
    times, dB_kam3, label3 = swpc.dB_kam_fromfile('dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt', fullname=False)
    assert(times.shape[0] == 31)

    assert(dB_kam1.shape[0] == 31)
    assert(dB_kam2.shape[0] == 31)


    times = np.array(times, dtype=int)
    dtimes = []
    for i in range(times.shape[0]):
        dtimes.append(datetime.datetime(times[i,0],times[i,1],times[i,2],times[i,3],times[i,4],times[i,5]))

    dB_kam_norm1 = np.sqrt(np.einsum('ij,ij->i', dB_kam1, dB_kam1))
    dB_kam_norm2 = np.sqrt(np.einsum('ij,ij->i', dB_kam2, dB_kam2))
    dB_kam_norm3 = np.sqrt(np.einsum('ij,ij->i', dB_kam3, dB_kam3))


    time_common, filenames, dB_SWMF_allcommon = swpc.commonTimes()
    dB_SWMF = dB_SWMF_allcommon[swpc.samp, :]
    dB_SWMF_norm = np.sqrt(np.einsum('ij,ij->i', dB_SWMF, dB_SWMF))


    plt.plot(dtimes, dB_kam_norm2-dB_kam_norm1, label='(d=1/4)-(d=1/2)', alpha=0.7)
    plt.plot(dtimes, dB_kam_norm3-dB_kam_norm2, label='(d=1/8)-(d=1/4)', alpha=0.7)
    plt.plot(dtimes, dB_kam_norm3-dB_SWMF_norm, label='(d=1/8)-SWMF', alpha=0.7)

    datetick('x')
    plt.xlabel('xlabel here')
    plt.title('dB, d comparison, range $\pm 32$')
    plt.legend()
    plt.show()

def plot_explore2():
    from hapiclient.plot.datetick import datetick
    import matplotlib.pyplot as plt
    import datetime

    times, dB_kam1, label1 = swpc.dB_kam_fromfile('dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.50000_.txt', fullname=False)
    assert(times.shape[0] == 31)
    times, dB_kam2, label2 = swpc.dB_kam_fromfile('dB_kam_tofile_-112.00_0016.00_-064.00_0064.00_-064.00_0064.00_0.50000_.txt', fullname=False)
    assert(times.shape[0] == 31)
    times, dB_kam3, label3 = swpc.dB_kam_fromfile('dB_kam_tofile_-224.00_0032.00_-128.00_0128.00_-128.00_0128.00_0.50000_.txt', fullname=False)
    assert(times.shape[0] == 31)

    assert(dB_kam1.shape[0] == 31)
    assert(dB_kam2.shape[0] == 31)


    times = np.array(times, dtype=int)
    dtimes = []
    for i in range(times.shape[0]):
        dtimes.append(datetime.datetime(times[i,0],times[i,1],times[i,2],times[i,3],times[i,4],times[i,5]))

    dB_kam_norm1 = np.sqrt(np.einsum('ij,ij->i', dB_kam1, dB_kam1))
    dB_kam_norm2 = np.sqrt(np.einsum('ij,ij->i', dB_kam2, dB_kam2))
    dB_kam_norm3 = np.sqrt(np.einsum('ij,ij->i', dB_kam3, dB_kam3))


    time_common, filenames, dB_SWMF_allcommon = swpc.commonTimes()
    dB_SWMF = dB_SWMF_allcommon[swpc.samp, :]
    dB_SWMF_norm = np.sqrt(np.einsum('ij,ij->i', dB_SWMF, dB_SWMF))


    plt.plot(dtimes, dB_kam_norm1-dB_SWMF_norm, label='($\pm 32$)-SWMF', alpha=0.7)
    plt.plot(dtimes, dB_kam_norm2-dB_SWMF_norm, label='($\pm 64$)-SWMF', alpha=0.7)
    plt.plot(dtimes, dB_kam_norm3-dB_SWMF_norm, label='($\pm 128$)-SWMF', alpha=0.7)

    datetick('x')
    plt.xlabel('xlabel here')
    plt.title('dB, range comparison, d=0.5')
    plt.legend()
    plt.show()


if plot:
    compsyst = 'MAG'
    swpc.plot_from_file(['dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt'], component='down', compsyst = compsyst)
    swpc.plot_from_file(['dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt'], component='east', compsyst = compsyst)
    swpc.plot_from_file(['dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt'], component='north', compsyst = compsyst)
    swpc.plot_from_file(['dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt'], component='norm', compsyst = compsyst)


    '''
    swpc.plot_from_file(['dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.50000_.txt',
                'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.25000_.txt',
                'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt'])

    swpc.plot_from_file(['dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.50000_.txt',
                'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.25000_.txt',
                'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt',
                'dB_kam_tofile_-112.00_0016.00_-064.00_0064.00_-064.00_0064.00_0.25000_.txt',
                'dB_kam_tofile_-112.00_0016.00_-064.00_0064.00_-064.00_0064.00_0.50000_.txt',
                'dB_kam_tofile_-224.00_0032.00_-128.00_0128.00_-128.00_0128.00_0.50000_.txt'])
    '''
    #plot_explore()

else:
    to = tm.time()

    time_common, filenames, dummy = swpc.commonTimes()
    ret = swpc.dB_kam_tofile(time_common, filenames, tag=None, xlims=(-112., 16.), ylims=(-64., 64.), zlims=(-64., 64.), d=0.25) #0.25

    tf = tm.time()


    print('time = ' + str((tf-to)/3600.) + ' hours') # time = 3.4527500342 hours
    g = open('done_time.txt', 'a')
    g.write('\n time = ' + str((tf-to)/3600.) + ' hours\n')
    g.close()
    print(ret)

    ###################### two runs
    to = tm.time()

    time_common, filenames, dummy = swpc.commonTimes()
    ret = swpc.dB_kam_tofile(time_common, filenames, tag=None, xlims=(-224., 32.), ylims=(-128., 128.), zlims=(-128., 128.), d=0.5)

    tf = tm.time()

    print('time = ' + str((tf-to)/3600.) + ' hours') # time = 3.4527500342 hours
    g = open('done_time.txt', 'a')
    g.write('\n time = ' + str((tf-to)/3600.) + ' hours\n')
    g.close()
    print(ret)

