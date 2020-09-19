import os
import sys
import numpy as np
import tempfile


from config import conf
import biot_savart_kameleon_interpolated_grid as bsk
import util
import cxtransform as cx
import magnetometers as mg



def getquadrants(run, time, location, fwrite=True): # loc in MAG sph

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

    regions = [q1, q2, q3, q4, q5, q6, q7, q8]

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

def main(run, location): # loc in MAG sph

    nf = 2
    para = False

    files = list(util.get_available_slices(run)[0])
    if nf is not None:
        files = files[:nf]

    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is not None and num_cores > len(files):
            num_cores = len(files)
        print('Parallel processing {0:d} file(s) using {1:d} cores'\
              .format(len(files), num_cores))
        values = Parallel(n_jobs=num_cores)(\
            delayed(getquadrants)(run, util.CDFfilename2time(run, filename), location, fwrite=False) for filename in files)

    else:
        values = []
        for filename in files:
            time = util.CDFfilename2time(run, filename)
            values.append(getquadrants(run, time, location, fwrite=False))

        values = np.array(values)

    assert(len(values.shape) == 3)
    outname = 'quadrants_values_%dx%dx%d.bin'%values.shape

    values.tofile(outname)

    print(values)
    print(values.shape)
    print('DONE')


if __name__=='__main__':
    run = 'DIPTSUR2'
    #time = (2019,9,2,6,30,0)
    location = mg.GetMagnetometerLocation('colaba', (2019,1,1,1,0,0), 'MAG', 'sph')
    main(run, location)
