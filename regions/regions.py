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

    def compute(region):
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

        print(np.max(np.abs(positive + negative - deltaB_loc)))
        assert(np.max(np.abs(positive + negative - deltaB_loc)) < 1e-9)

        # index k runs from: 
                # k=0 -> 'positive'
                # k=1 -> 'negative'
                # k=2 -> 'deltaB_loc' = 'positive' + 'negative'
        # index l runs from: 
                # l=0 -> 'north'
                # l=1 -> 'east'
                # l=2 -> 'down'
        return [positive, negative, deltaB_loc] # indexed by above (k,l)

    toret = []
    for i in range(len(regions)):
        ret = compute(regions[i])
        toret.append(ret)

        if fwrite:
            if os.path.exists('/home/gary/'):
                f = open('/home/gary/temp/' + run + 'regs.txt','a')
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

    # index j runs from:
            # j = 0 -> 1st region
            # ...
            # j = n-1 -> nth region
    # index k runs from: 
            # k=0 -> 'positive'
            # k=1 -> 'negative'
            # k=2 -> 'deltaB_loc' = 'positive' + 'negative'
    # index l runs from: 
            # l=0 -> 'north'
            # l=1 -> 'east'
            # l=2 -> 'down'
    return np.array(toret) # indexed by above (j,k,l)


def plot(run, pkl, comp, show=False, tag='', totxt=False):
    if isinstance(comp, int):
        icomp = comp
        comp = str(comp)
    elif comp=='north':
        icomp = 0
    elif comp=='east':
        icomp = 1
    elif comp=='down':
        icomp = 2
    else:
        raise ValueError ("component must be 'north', 'east', 'down', or integer")

    pkl = conf[run+'_derived'] + 'regions/' + pkl

    with open(pkl, 'rb') as handle:
        result = pickle.load(handle)

    #result['shape_README']
    location = result['location']
    regions  = result['regions']
    deltaBs = result['deltaBs']
    times = result['times']

    assert(deltaBs.shape[0] == times.shape[0])

    timesfull = util.get_available_slices(run)[1]
    if not np.all(timesfull == times):
        print('WARNING: appears to be time missmatch. Might not include all files of run')

    #if nf is not None:
    #    times = times[:nf, :]
    #assert(times.shape[0] == shape[0])

    import datetime
    dtimes = []
    for i in range(times.shape[0]):
        dtimes.append(datetime.datetime(times[i,0],times[i,1],times[i,2],times[i,3],times[i,4],times[i,5]))

    types = ('positive', 'negative', 'deltaB_loc') #!!!!!!!!!!
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
            plt.plot(dtimes, deltaBs[:, i, j, icomp], label=types[j])
            if totxt:
                txt.write('\nlabel=%s\n'%(types[j]))
                np.savetxt(txt, np.column_stack([ times, deltaBs[:, i, j, icomp] ]), fmt='%.5f')
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


def signedintegrate_timeseries(run, location, regions='octants', tag=''): # location in MAG sph

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
        deltaBs = Parallel(n_jobs=num_cores)(\
            delayed(signedintegrate)(run, time, location, regions=regions) for time in list(times))

    else:
        deltaBs = []
        for time in list(times):
            deltaBs.append(signedintegrate(run, time, location, regions=regions))

    # index i runs from:
            # i = 0 -> 1st time there's a datafile in cdflist.txt
            # ...
            # i = N-1 -> Nth time there's a datafile in cdflist.txt
    # index j runs from:
            # j = 0 -> 1st region
            # ...
            # j = n-1 -> nth region
    # index k runs from: 
            # k=0 -> 'positive'
            # k=1 -> 'negative'
            # k=2 -> 'deltaB_loc' = 'positive' + 'negative'
    # index l runs from: 
            # l=0 -> 'north'
            # l=1 -> 'east'
            # l=2 -> 'down'
    deltaBs = np.array(deltaBs) # indexed by above (i,j,k,l)
    assert(len(deltaBs.shape) == 4)
    assert(deltaBs.shape[0] == times.shape[0])

    README_string = \
    ('"deltaBs" indexed (i,j,k,l) as follows:\n'
     ' index i runs from:\n'
     '       i = 0 -> 1st time there is a datafile in cdflist.txt\n'
     '       ...\n'
     '       i = N-1 -> Nth time there is a datafile in cdflist.txt\n'
     ' index j runs from:\n'
     '       j = 0 -> 1st region\n'
     '       ...\n'
     '       j = n-1 -> nth region\n'
     ' index k runs from:\n'
     '       k=0 -> positive_contribution\n'
     '       k=1 -> negative_contribution\n'
     '       k=2 -> full_deltaB = positive_contribution + negative_contribution\n'
     ' index l runs from:\n'
     '       l=0 -> north component\n'
     '       l=1 -> east component\n'
     '       l=2 -> down component\n'
     '\n'
     '"times" indexed (i,j) as follows:\n'
     ' index i runs from:\n'
     '       i = 0 -> 1st time there is a datafile in cdflist.txt\n'
     '       ...\n'
     '       i = N-1 -> Nth time there is a datafile in cdflist.txt\n'
     ' index j runs from:\n'
     '       j = 0 -> year\n'
     '       j = 1 -> month\n'
     '       j = 2 -> day\n'
     '       j = 3 -> hours\n'
     '       j = 4 -> minutes\n'
     '       j = 5 -> seconds\n'
     '       j = 6 -> miliseconds\n')

    result =   {'README' : README_string,
                'location' : location,
                'regions' : regions,
                'deltaBs' : deltaBs,
                'times' : times
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

    print(deltaBs.shape)
    print('DONE')

    return pkl


def main():
    run = 'IMP10_RUN_SAMPLE'
    time = (2019,9,2,7,0,0)
    location = mg.GetMagnetometerLocation('colaba', (2019,1,1,1,0,0), 'MAG', 'sph')


    pm = 31.875
    reg =  {'xlims': (-pm, pm),
            'ylims': (-pm, pm),
            'zlims': (-pm, pm),
            'd': 0.25
            }
    signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=0.)
    signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=1.475)
    signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=1.525)
    signedintegrate(run, time, location, regions=(reg,), fwrite=True, rmin=1.55)

    if False:
        location = mg.GetMagnetometerLocation('colaba', (2019,1,1,1,0,0), 'MAG', 'sph')
        signedintegrate_timeseries(run, location, regions='octants', tag='octants') 
    else:
        #plot(run, 'mlat_11.059_mlon_146.897_nf_240-octants.pkl', 'north', tag='north', totxt=True)
        #plot(run, 'mlat_11.017_mlon_147.323_nf_698-octants.pkl', 'north', tag='north', totxt=True)
        print('else')

if __name__=='__main__':
    main()