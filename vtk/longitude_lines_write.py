# longitude_lines_write

import sys
import os
import numpy as np

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from util import maketag
import cxtransform as cx

# run parameters
#Nlong = 5
#Nb = 6
sign = -1  # changes sign of magnetic field used to trace the field lines
debug = True

'''

def write_line_vtk(Event, lon_array=[0., 10., -10., 20., -20.]):

    if isinstance(Event[0],list):
        time = Event[1]
    else:
        time = Event[0:5]

    tag = maketag(time)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(ret['run_path_derived'])
    for lon in lon_array:
        out_fname = conf["run_path_derived"] + subdir + 'GEO_longitude_line_%.2f.vtk' %(lon,)
        lat = np.linspace(-90,90,100)
        a = np.column_stack([np.ones(100, ), lat, lon*np.ones(100, )])
        sol = cx.GEOtoGSM(a, time, 'sph', 'car')

        writevtk(out_fname, sol, None, 'line', None, title='Title', ftype='BINARY', grid='POLYDATA')

'''

def write_line(sol, out_fname):
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

if True:
    time = [2000, 1, 1, 0, 0, 0]
    #tag = maketag(time)
    #subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    #if not os.path.exists(conf["run_path_derived"] + subdir):
    #    os.mkdir(ret['run_path_derived'])

    for lon in [0.]:
        out_fname = '/tmp/GEO_other_longitude_line_%.2f.vtk' %(lon,)
        lat = np.linspace(-90,90,100)
        a = np.column_stack([np.ones(100, ), lat, lon*np.ones(100, )])
        sol = cx.GEOtoGSM(a, time, 'sph', 'car')

        write_line(sol, out_fname)
        #writevtk(out_fname, sol, None, 'line', None, title='Title', ftype='ASCII', grid='POLYDATA')

