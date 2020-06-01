# longitude_lines_write

import sys
import os
import numpy as np

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import pos_sun as ps
import cut_plane

from scipy.integrate import odeint
import _CCMC as ccmc

# run parameters
#Nlong = 5
#Nb = 6
sign = -1  # changes sign of magnetic field used to trace the field lines
debug = True

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

def writevtk(Event, lon_array=[0., 10., -10., 20., -20.]):
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    time = Event[0:6]
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d' % tuple(time)

    #solns_restr = Compute(Event, Nb)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(ret['run_path_derived'])
    for lon in lon_array:
        out_fname = conf["run_path_derived"] + subdir + 'GEO_longitude_line_%.2f.vtk' %(lon,)
        lat = np.linspace(-90,90,100)
        a = np.column_stack([np.ones(100, ), lat, lon*np.ones(100, )])
        sol = ps.GEOtoGSM(a, time, 'sph', 'car')
        if debug:
            print('Writing ' + out_fname)
        f = open(out_fname, 'w')
        f.write('# vtk DataFile Version 3.0\n')
        f.write('A dataset with one polyline and no attributes\n')
        f.write('ASCII\n')
        f.write('\n')
        f.write('DATASET POLYDATA\n')
        f.write('POINTS ' + str(sol.shape[0]) + ' float\n')
        for k in range(sol.shape[0]):
            f.write('%e %e %e\n' % (sol[k, 0], sol[k, 1], sol[k, 2]))
        f.write('LINES ' + '1' + ' ' + str(sol.shape[0] + 1) + '\n' )
        f.write(str(sol.shape[0]) + '\n')
        for k in range(sol.shape[0]):
            f.write(str(k) + '\n')
        f.close()
        if debug:
            print('Wrote ' + out_fname)
