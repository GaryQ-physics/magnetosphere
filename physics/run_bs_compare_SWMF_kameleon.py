import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import bs_compare_SWMF_kameleon as bscomp
#import cxtransform as cx


plot = True
run = 'SCARR5'


if plot:
    #compsyst = 'MAG'
    bscomp.plot_from_file(run, ['dB_SWMF_tofile.txt',
     'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.50000_.txt', 
     'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.12500_.txt'])

else:
    bscomp.compute(run, tag=None, xlims=(-112., 16.), ylims=(-64., 64.), zlims=(-64., 64.), d=0.25)

