import os
import sys
import numpy as np
import tempfile


from config import conf
import biot_savart_kameleon_interpolated_grid as bsk
import util
import cxtransform as cx
import magnetometers as mg


def getoctants():
    q1 = {'xlims': (0., 32.),
          'ylims': (0., 32.),
          'zlims': (0., 32.),
          'd': 0.25

            }

    q2 = {'xlims': (0., 32.),
          'ylims': (0., 32.),
          'zlims': (-32., 0.),
          'd': 0.25

            }

    q3 = {'xlims': (0., 32.),
          'ylims': (-32., 0),
          'zlims': (0., 32.),
          'd': 0.25

            }

    q4 = {'xlims': (-32., 0.),
          'ylims': (0., 32.),
          'zlims': (0., 32.),
          'd': 0.25

            }

    q5 = {'xlims': (0., 32.),
          'ylims': (-32., 0.),
          'zlims': (-32., 0.),
          'd': 0.25

            }

    q6 = {'xlims': (-32., 0.),
          'ylims': (0., 32.),
          'zlims': (-32., 0.),
          'd': 0.25

            }

    q7 = {'xlims': (-32., 0.),
          'ylims': (-32., 0.),
          'zlims': (0., 32.),
          'd': 0.25

            }

    q8 = {'xlims': (-32., 0.),
          'ylims': (-32., 0.),
          'zlims': (-32., 0.),
          'd': 0.25

            }

    return [q1, q2, q3, q4, q5, q6, q7, q8]


def signedintegrate(run, time, location, fwrite=True): # loc in MAG sph

    regions = getoctants()

    def bla(region):
        xlims = region['xlims']
        ylims = region['ylims']
        zlims = region['zlims']
        d = region['d']

        dB, G, Ntup = bsk.integrate(run, time, location[1], location[2], para=False,
            xlims=xlims, ylims=ylims, zlims=zlims, d=d, returnAll=True)
        deltaB = np.sum(dB, axis=0)
        deltaB_loc = bsk.toMAGLocalComponents(time, location[1], location[2], deltaB)

        dB_loc = bsk.toMAGLocalComponents(time, location[1], location[2], dB)

        positive = np.empty(3)
        negative = np.empty(3)
        for comp in range(3):
            tr = dB_loc[:, comp]>0
            positive_contrs = dB_loc[:,comp][tr]
            negative_contrs = dB_loc[:,comp][np.logical_not(tr)]
            positive[comp] =  np.sum(positive_contrs)
            negative[comp] =  np.sum(negative_contrs)

        #assert(np.max(np.abs(positive + negative - deltaB_loc))<1e-6)
        print(np.max(np.abs(positive + negative - deltaB_loc)))
        return [positive, negative, deltaB_loc]

    toret = []

    for region in regions:

        ret = bla(region)

        if fwrite:
            f = open('/home/gary/temp/quads.txt','a')
            #f = open('/media/solar-backup/tmp/quads.txt','a')

            f.write('\n\n')
            f.write('time = ' + str(time))
            f.write('\nmlat %f, mlon %f'%(location[1], location[2]))
            f.write('\noctant(GSM):\n' + str(region))
            f.write('\nnet positive contributions = '+str(ret[0])+'  (north, east, down)')
            f.write('\nnet negative contributions = '+str(ret[1])+'  (north, east, down)')
            f.write('\ndeltaB_loc/nT = ' + str(ret[2]))
            f.write('\n\n')

            f.close()

        toret.append(ret)

    return np.array(toret)


def compidx(comp):
    if comp=='north':
        return 0
    if comp=='east':
        return 1
    if comp=='down':
        return 2


def plot(filename, shape, comp='north'):

    times = util.get_available_slices(run)[1]
    assert(times.shape[0] == shape[0])

    import datetime
    dtimes = []
    for i in range(times.shape[0]):
        dtimes.append(datetime.datetime(times[i,0],times[i,1],times[i,2],times[i,3],times[i,4],times[i,5]))

    values = np.fromfile('quadrants_values_%dx%dx%dx%d.bin'%shape).reshape(shape)

    regions = getoctants()
    types = ['positive', 'negative', 'deltaB_loc'] #!!!!!!!!!!

    for i in range(len(regions)):
        title = 'dB_' + comp + str(regions[i])

        from hapiclient.plot.datetick import datetick

        import matplotlib.pyplot as plt
        import datetime

        for j in range(3):
            #print('plot add line   dtimes vs values[:, i, %d, compidx(comp)]'%(j))
            #print('with legend %s'%(types[j]))
            #print(dtimes)
            #print(values[:, i, j, compidx(comp)])
            plt.plot(dtimes, values[:, i, j, compidx(comp)], label=types[j])

        print('export png for %d with title %s'%(i,title))


        datetick('x')
        # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
        #plt.gcf().autofmt_xdate()

        plt.xlabel('time')
        plt.ylabel('(in nanoTesla)')

        plt.title(title)
        plt.legend()

        #plt.show()
        plt.savefig('quadrants_%d.png'%(i))
        plt.clf()


def main(run, location): # loc in MAG sph

    nf = None
    para = True

    files = list(util.get_available_slices(run)[0])
    times = util.get_available_slices(run)[1]

    if nf is not None:
        files = files[:nf]
        times = times[:nf, :]

    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is not None and num_cores > len(files):
            num_cores = len(files)
        print('Parallel processing {0:d} file(s) using {1:d} cores'\
              .format(len(files), num_cores))
        values = Parallel(n_jobs=num_cores)(\
            #delayed(signedintegrate)(run, util.CDFfilename2time(run, filename), location, fwrite=False) for filename in files)
            delayed(signedintegrate)(run, time, location, fwrite=False) for time in list(times))

    else:
        values = []
        #for filename in files:
        #    time = util.CDFfilename2time(run, filename)
        #    values.append(signedintegrate(run, time, location, fwrite=False))
        for time in list(times):
            values.append(signedintegrate(run, time, location, fwrite=False))


    values = np.array(values)
    print(values)
    assert(len(values.shape) == 4)
    direct = conf[run +'_derived'] + 'regions/'
    if not os.path.exists(direct):
        os.makedirs(direct)
    outname = direct + 'quadrants_values_%dx%dx%dx%d.bin'%values.shape

    util.safenumpy_tofile(values, outname)

    #print(values)
    print(values.shape)
    print('DONE')

    return outname


if __name__=='__main__':
    run = 'DIPTSUR2'
    #time = (2019,9,2,6,30,0)
    location = mg.GetMagnetometerLocation('colaba', (2019,1,1,1,0,0), 'MAG', 'sph')
    #main(run, location)
    plot('quadrants_values_698x8x3x3.bin', (698,8,3,3), comp='north')
