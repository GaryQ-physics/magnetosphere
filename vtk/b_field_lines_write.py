# b_field_lines_write

import sys
import os
import numpy as np

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import cutplane as cut
from util import time2filename, maketag
import cxtransform as cx




from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import odeint


def Compute(Event, Nb, debug=False):
    Rfudge = 1.01
    if isinstance(Event[0],list):
        mag = np.array(Event[0])
        time = Event[1]
    else:
        time = Event[0:5]
        mag = np.array([Rfudge, Event[5], Event[6]])

    filename = time2filename(time)

    R = mag[0]
    MLAT = mag[1]
    MLON = mag[2]

    eps = 3.
    #IC = []
    ICmag = []
    for i in range(Nb+1):
        if i==0:
            delta = 0.
        elif i>Nb/2:
            #delta=(i-Nb/2)*eps
            delta = (-i)*eps
        else:
            #delta=(i-Nb/2)*eps
            delta = (-i)*eps
        #IC.append(ps.MAGtoGSM([R, MLAT-delta, MLON], time[0:6], 'sph', 'car'))
        ICmag.append(np.array([R, MLAT-delta, MLON]))

    if debug:
        print(ICmag)

    solns_restr = cut.fieldlines(time, ICmag, debug=debug, fieldvar='b')
    if debug:
        print(len(solns_restr))
        print(solns_restr)

    return solns_restr


def writevtk(Event, Nb=6, debug=True):
    solns_restr = Compute(Event, Nb)

    if isinstance(Event[0],list):
        mag = np.array(Event[0])
        time = Event[1]
    else:
        time = Event[0:5]
        mag = np.array([1., Event[5], Event[6]])

    tag = maketag(time)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)
    for i in range(Nb+1):
        from_list = solns_restr[i]
        sol = np.array(from_list)
        v = cx.GSMtoMAG([sol[0, 0], sol[0, 1], sol[0, 2]], time[0:6], 'car', 'sph')
        if i==0:
            out_fname = conf["run_path_derived"] + subdir + 'b' +'_field_line_event_%.2f_%.2f.vtk' %(mag[1], mag[2]) #(MLAT, MLON)
        else:
            out_fname = conf["run_path_derived"] + subdir + 'b' +'_field_line_%.2f_%.2f.vtk' %(v[1], v[2]) #(MLAT, MLON)

        if debug: 
            print('Writing ' + out_fname)
            if i==0:
                quest = '%.2f_%.2f' %(v[1], v[2]) == '%.2f_%.2f' %(mag[1], mag[2]) # should be True
                print('should be true: ' + str(quest))
        f = open(out_fname, 'w')
        f.write('# vtk DataFile Version 3.0\n')
        f.write('A dataset with one polyline and no attributes\n')
        f.write('ASCII\n')
        f.write('\n')
        f.write('DATASET POLYDATA\n')
        f.write('POINTS ' + str(sol.shape[0]) + ' float\n')
        for k in range(sol.shape[0]):
            f.write('%e %e %e\n'%(sol[k, 0], sol[k, 1], sol[k, 2]))
        f.write('LINES ' + '1' + ' ' + str(sol.shape[0] + 1) + '\n' )
        f.write(str(sol.shape[0]) + '\n')
        for k in range(sol.shape[0]):
            f.write(str(k) + '\n')
        f.close()
        if debug:
            print('Wrote ' + out_fname)
