import os
import sys
import numpy as np
import tempfile


from config import conf
import biot_savart_kameleon_interpolated_grid as bsk
import util
import cxtransform as cx
import magnetometers as mg



def main(run, time, location): # loc in MAG sph

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

    for region in regions:

        ret = bla(region)

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



if __name__=='__main__':
    run = 'DIPTSUR2'
    time = (2019,9,2,6,30,0)
    location = mg.GetMagnetometerLocation('colaba', time, 'MAG', 'sph')
    main(run, time, location)
