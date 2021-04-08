import numpy as np
from numba import njit
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
from read_swmf_files2 import read_all


@njit
def _jit_probe():
    # needed ?
    pass


def slice_probe(run, time, obs_point, cache=None):
    if cache is None:
        cache = read_all(util.time2CDFfilename(run,time)[:-8])

    ###########
    assert(point[0] == NeededArray[0,9477, 3,1,4])
    assert(point[1] == NeededArray[1,9477, 3,1,4])
    assert(point[2] == NeededArray[2,9477, 3,1,4])

    b1x_nat = NeededArray[_b1x,9477, 3,1,4]
    b1y_nat = NeededArray[_b1y,9477, 3,1,4]
    b1z_nat = NeededArray[_b1z,9477, 3,1,4]
    ###########


def stitch_probe(run, time, obs_point):
    pass
