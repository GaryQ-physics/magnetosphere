import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import bs_compare_SWMF_kameleon as bscomp
#import cxtransform as cx


plot = True

#run = 'SWPC'
#station = 'YKC'

run = 'SCARR5'
station = [63.362, 0.] # aproximatly step of 175./174

if plot:
    #compsyst = 'MAG'
    bscomp.plot_from_file(run, station, ['dB_SWMF_tofile.txt',
     'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.50000_.txt',
     #'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt',
     #'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_1.00000_.txt',
        ])

else:
    bscomp.compute(run, station, tag=None, xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.125, skip_SWMF=False)

