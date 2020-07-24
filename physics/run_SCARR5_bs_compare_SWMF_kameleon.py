import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import SCARR5_bs_compare_SWMF_kameleon as scarr5
#import cxtransform as cx
import time as tm


plot = True

if plot:
    compsyst = 'MAG'
    scarr5.plot_from_file(['dB_kam_tofile_-048.00_0016.00_-032.00_0032.00_-032.00_0032.00_0.50000_.txt'])

else:
    to = tm.time()

    time_common, filenames, dummy = scarr5.commonTimes()
    ret = scarr5.dB_kam_tofile(time_common, filenames, tag=None, xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.5)

    tf = tm.time()


    print('time = ' + str((tf-to)/3600.) + ' hours') # time = 3.4527500342 hours
    g = open('done_time.txt', 'a')
    g.write('\n time = ' + str((tf-to)/3600.) + ' hours\n')
    g.close()
    print(ret)

