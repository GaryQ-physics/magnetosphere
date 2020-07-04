import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import cxtransform as cx



def writevtk(Event, Nt=100, Np=100):
    #Event = [year, month, day, hours, minutes, seconds, MLONdeg, MLATdeg]
    time = Event[0:6]
    tag = '_%04d:%02d:%02dT%02d:%02d:%02d' % tuple(time)
    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])

    fname = conf["run_path_derived"] + subdir + 'earth' + tag +'.vtk'

    R = 1.
    theta = np.linspace(0., np.pi, Nt)
    phi = np.linspace(0., 2.*np.pi, Np)

    B1, B2 = np.meshgrid(phi, theta)
    B1= B1.flatten(order='C')
    B2= B2.flatten(order='C')
    B = np.column_stack((B1, B2))  # B[i, 0]=B1[i], B[i, 1]=B2[i]

    normPhi = np.linspace(0., 1., Np)
    normTheta = np.flipud(np.linspace(0., 1., Nt))
    u, v = np.meshgrid(normPhi, normTheta)
    u = u.flatten(order='C')
    v = v.flatten(order='C')
    UV = np.column_stack((u, v))

    PI = np.pi*np.ones((B1.size, ))
    x = R*np.cos(B1+PI)*np.sin(B2)
    y = R*np.sin(B1+PI)*np.sin(B2)
    z = R*np.cos(B2)
    XYZ = np.column_stack((x, y, z))

    XYZr = cx.GEOtoGSM(XYZ, time, 'car', 'car')

    print("Writing " + fname)
    f = open(fname,'w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('Structured Grid for rotated sphere\n')
    f.write('ASCII\n')
    f.write('DATASET STRUCTURED_GRID\n')
    f.write('DIMENSIONS ' + str(Nt) + ' ' + str(Np) + ' ' + str(1) + '\n' )
    f.write('POINTS '+str(Nt*Np)+' float\n')
    np.savetxt(f, XYZr)
    f.write('\n')

    f.write('POINT_DATA ' + str(Nt*Np) + '\n')
    f.write('TEXTURE_COORDINATES TextureCoordinates 2 float\n') # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping
    np.savetxt(f, UV)

    f.close()
    print("Wrote " + fname)
