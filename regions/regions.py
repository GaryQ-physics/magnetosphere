import os
import sys
import numpy as np
import tempfile
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import biot_savart_kameleon_interpolated_grid as bsk
import util
import cxtransform as cx
import magnetometers as mg

#https://stackoverflow.com/questions/16779497/how-to-set-memory-limit-for-thread-or-process-in-python
import resource
soft, hard = 32*10**30, 32*10**30
resource.setrlimit(resource.RLIMIT_AS,(soft, hard))

def getfull(pm=16.):
    q = {'xlims': (-pm, pm),
         'ylims': (-pm, pm),
         'zlims': (-pm, pm),
         'd': 0.25
            }

    return (q,)

def getoctants():
    q1 = {'xlims': (0., 32.),
          'ylims': (0., 32.),
          'zlims': (0., 32.),
          'd': 0.25

            }

    q2 = {'xlims': (0., 32.),
          'ylims': (0., 32.),
          'zlims': (-32., -0.25),
          'd': 0.25

            }

    q3 = {'xlims': (0., 32.),
          'ylims': (-32., -0.25),
          'zlims': (0., 32.),
          'd': 0.25

            }

    q4 = {'xlims': (-32., -0.25),
          'ylims': (0., 32.),
          'zlims': (0., 32.),
          'd': 0.25

            }

    q5 = {'xlims': (0., 32.),
          'ylims': (-32., -0.25),
          'zlims': (-32., -0.25),
          'd': 0.25

            }

    q6 = {'xlims': (-32., -0.25),
          'ylims': (0., 32.),
          'zlims': (-32., -0.25),
          'd': 0.25

            }

    q7 = {'xlims': (-32., -0.25),
          'ylims': (-32., -0.25),
          'zlims': (0., 32.),
          'd': 0.25

            }

    q8 = {'xlims': (-32., -0.25),
          'ylims': (-32., -0.25),
          'zlims': (-32., -0.25),
          'd': 0.25

            }

    return (q1, q2, q3, q4, q5, q6, q7, q8)


def signedintegrate(run, time, location, regions='octants', fwrite=False, rmin=None): # loc in MAG sph

    if regions == 'octants':
        regions = getoctants()
    elif regions == 'full':
        regions = getfull()
    elif isinstance(regions, str):
        raise ValueError ("regions must be 'octants', 'full', or custom tuple of dictionaries")

    def bla(region):
        xlims = region['xlims']
        ylims = region['ylims']
        zlims = region['zlims']
        d = region['d']

        dB, G, Ntup = bsk.integrate(run, time, location[1], location[2], para=False,
            xlims=xlims, ylims=ylims, zlims=zlims, d=d, returnAll=True, rmin=rmin)
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
        assert(np.max(np.abs(positive + negative - deltaB_loc)) < 1e-9)
        return [positive, negative, deltaB_loc]

    toret = []

    for i in range(len(regions)):
        ret = bla(regions[i])
        toret.append(ret)

        if fwrite:
            if os.path.exists('/home/gary/'):
                f = open('/home/gary/temp/regs.txt','a')
            else:
                f = open('/tmp/regs.txt','a')

            f.write('\n\n')
            f.write('run = %s\n'%(run))
            f.write('time = ' + str(time))
            f.write('\nmlat %f, mlon %f'%(location[1], location[2]))
            f.write('\noctant(GSM):with rmin '+str(rmin)+':\n' + str(regions[i]))
            f.write('\nnet positive contributions = '+str(ret[0])+'  (north, east, down)')
            f.write('\nnet negative contributions = '+str(ret[1])+'  (north, east, down)')
            f.write('\ndeltaB_loc/nT = ' + str(ret[2]))
            f.write('\n\n')

            f.close()


    return np.array(toret)


def compidx(comp):
    if comp=='north':
        return 0
    if comp=='east':
        return 1
    if comp=='down':
        return 2
    raise ValueError ("component must be 'north', 'east', or 'down'")


def plot(run, pkl, comp, show=False, tag='', totxt=False):
    pkl = conf[run+'_derived'] + 'regions/' + pkl

    with open(pkl, 'rb') as handle:
        result = pickle.load(handle)

    nf = result['nf']
    shape = result['shape']
    location = result['location']
    regions  = result['regions']
    values = result['values']
    assert(values.shape == shape)
    #values = values.reshape(shape)

    times = util.get_available_slices(run)[1]
    if nf is not None:
        times = times[:nf, :]
    assert(times.shape[0] == shape[0])

    import datetime
    dtimes = []
    for i in range(times.shape[0]):
        dtimes.append(datetime.datetime(times[i,0],times[i,1],times[i,2],times[i,3],times[i,4],times[i,5]))

    types = ['positive', 'negative', 'deltaB_loc'] #!!!!!!!!!!
    for i in range(len(regions)):
        title = 'dB_%s_'%(comp) + 'mlat_{0:.3f}_mlon_{1:.3f}'.format(location[1], location[2]) \
                +'\n' +str(regions[i])
        outname = conf[run+'_derived'] + 'regions/%s_region_%d.png'%(tag,i)

        from hapiclient.plot.datetick import datetick

        import matplotlib.pyplot as plt
        import datetime

        if totxt:
            print('writing txt ' + outname + '.txt')
            txt = open(outname + '.txt', 'w')
            txt.write('title is:\n')
            txt.write(title +'\n')
            txt.write('year month day hour minute second milisecond value(label)\n')
        for j in range(3):
            #print('plot add line   dtimes vs values[:, i, %d, compidx(comp)]'%(j))
            #print('with legend %s'%(types[j]))
            #print(dtimes)
            #print(values[:, i, j, compidx(comp)])
            plt.plot(dtimes, values[:, i, j, compidx(comp)], label=types[j])
            if totxt:
                txt.write('\nlabel=%s\n'%(types[j]))
                np.savetxt(txt, np.column_stack([ times, values[:, i, j, compidx(comp)] ]), fmt='%.5f')
        print('export png for %d with title %s'%(i,title))


        datetick('x')
        # https://stackoverflow.com/questions/1574088/plotting-time-in-python-with-matplotlib
        #plt.gcf().autofmt_xdate()

        plt.xlabel('time')
        plt.ylabel('(in nanoTesla)')

        plt.title(title)
        plt.legend()

        if show: plt.show()
        plt.savefig(outname)
        plt.clf()
        if totxt: txt.close()


def main(run, location, regions='octants', tag=''): # loc in MAG sph

    if regions == 'octants':
        regions = getoctants()
    elif regions == 'full':
        regions = getfull()
    elif isinstance(regions, str):
        raise ValueError ("regions must be 'octants', 'full', or custom tuple of dictionaries")

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
            delayed(signedintegrate)(run, time, location, regions=regions) for time in list(times))

    else:
        values = []
        #for filename in files:
        #    time = util.CDFfilename2time(run, filename)
        #    values.append(signedintegrate(run, time, location, fwrite=False))
        for time in list(times):
            values.append(signedintegrate(run, time, location, regions=regions))

    values = np.array(values)
    assert(len(values.shape) == 4)

    result =   {'nf' : nf,
                'shape' : values.shape,
                'location' : location,
                'regions' : regions,
                'values' : values
                }

    direct = conf[run +'_derived'] + 'regions/'
    if not os.path.exists(direct):
        os.makedirs(direct)

    if nf is None:
        nfstr = str(len(files))
    else:
        nfstr = str(nf)
    pkl = direct + 'mlat_{0:.3f}_mlon_{1:.3f}_nf_{2:s}-{3:s}.pkl' \
                    .format(location[1], location[2], nfstr, tag)

    print('writing to ' + pkl)
    util.safeprep_fileout(pkl)
    with open(pkl, 'wb') as handle:
        pickle.dump(result, handle)

    print(values.shape)
    print('DONE')

    return pkl


if __name__=='__main__':
    run = 'IMP10_RUN_SAMPLE'
    time = (2019,9,2,7,0,0)
    location = mg.GetMagnetometerLocation('colaba', (2019,1,1,1,0,0), 'MAG', 'sph')

    pm = 32.
    reg =  {'xlims': (-pm, pm),
            'ylims': (-pm, pm),
            'zlims': (-pm, pm),
            'd': 0.25
            }
    print(signedintegrate(run, time, location, regions=(reg,), fwrite=False, rmin=0.))

    if False:

        pm = 32.
        reg =  {'xlims': (-pm, pm),
                'ylims': (-pm, pm),
                'zlims': (-pm, pm),
                'd': 0.25
                }
        signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=1.)

        pm = 32.
        reg =  {'xlims': (-pm, pm),
                'ylims': (-pm, pm),
                'zlims': (-pm, pm),
                'd': 0.25
                }
        signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=1.5)

        pm = 32.
        reg =  {'xlims': (-pm, pm),
                'ylims': (-pm, pm),
                'zlims': (-pm, pm),
                'd': 0.25
                }
        signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=1.75)

        pm = 31.875
        reg =  {'xlims': (-pm, pm),
                'ylims': (-pm, pm),
                'zlims': (-pm, pm),
                'd': 0.25
                }
        signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=0.)

        pm = 31.875
        reg =  {'xlims': (-pm, pm),
                'ylims': (-pm, pm),
                'zlims': (-pm, pm),
                'd': 0.25
                }
        signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=1.)

        pm = 31.875
        reg =  {'xlims': (-pm, pm),
                'ylims': (-pm, pm),
                'zlims': (-pm, pm),
                'd': 0.25
                }
        signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=1.5)

        pm = 31.875
        reg =  {'xlims': (-pm, pm),
                'ylims': (-pm, pm),
                'zlims': (-pm, pm),
                'd': 0.25
                }
        signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=1.75)

        #time = (2019,9,2,6,30,0)
        signedintegrate(run, time, location, regions='octants', fwrite=True)
        #signedintegrate(run, time, location, regions='full', fwrite=True)


        pm = 16.
        reg =  {'xlims': (-pm, pm),
                'ylims': (-pm, pm),
                'zlims': (-pm, pm),
                'd': 0.25
                }

        signedintegrate(run, time, location, regions=(reg,), fwrite=True)


        pm = 32.
        reg =  {'xlims': (-pm, pm),
                'ylims': (-pm, pm),
                'zlims': (-pm, pm),
                'd': 0.25
                }

        signedintegrate(run, time, location, regions=(reg,), fwrite=True)

        #time = (2019,9,2,6,30,0)
        location = mg.GetMagnetometerLocation('colaba', (2019,1,1,1,0,0), 'MAG', 'sph')
        main(run, location, regions='octants', tag='octants') 
    else:
        #plot(run, 'mlat_11.059_mlon_146.897_nf_240-octants.pkl', 'north', tag='north', totxt=True)
        #plot(run, 'mlat_11.017_mlon_147.323_nf_698-octants.pkl', 'north', tag='north', totxt=True)
        print('else')

