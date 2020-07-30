import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import bs_compare_SWMF_kameleon as bscomp
#import cxtransform as cx


plot = True

run = 'SWPC'
station = 'FUR'

#run = 'SCARR5'
#station = [63.362, 0.] # aproximatly step of 175./174

if plot:

    #compsyst = 'MAG'
    bscomp.plot_from_file(run, station, ['dB_SWMF_tofile.txt',
     #'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_1.00000_.txt',
     'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.50000_.txt',
     #'dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.06250_.txt',
     #'dB_kam_tofile_-112.00_0016.00_-064.00_0064.00_-064.00_0064.00_0.12500_.txt',
     ])
    
    #bscomp.plot_from_file(run, station, ['/home/gary/magnetosphere/data/SWPC_SWMF_052811_2-derived/YKC/dB_SWMF_tofile.txt', '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2-derived/FRN/dB_SWMF_tofile.txt'], fullname=True)
    
else:
    bscomp.compute(run, station, tag=None, xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.5, skip_SWMF=False)
    #bscomp.compute(run, station, tag=None, xlims=(-112., 16.), ylims=(-64., 64.), zlims=(-64., 64.), d=0.125, skip_SWMF=True)

