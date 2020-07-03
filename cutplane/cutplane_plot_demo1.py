import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import cutplane_plot as cp
import cxtransform as cx

time = [2003, 11, 20, 7, 0, 0, 0]

mag = cx.GSMtoMAG([-1, 0, 0.7], time, 'car', 'sph')
#Mdipole = cx.MAGtoGSM([0., 0., 1.], time[0:6], 'car', 'car')
#print(Mdipole)


cp.plot(time, 'p', 'xz', logz=True, dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])
#cp.plot(time, 'p', mag, field_lines=field_lines, dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])
#cp.plot(time, 'p', 'xz', field_lines=mag, dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])
#cp.plot(time, 'p', [[1, 0, 0], [0, 0, 1]], field_lines=mag, dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])

#mag = cx.GSMtoMAG([-1, 0, 0.7], time, 'car', 'sph')
#cp.plot(time, 'p', mag, dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])

if True:
    dtheta = 0.5 # takes a while with 0.1 with probe method
    import numpy as np
    angles = (np.pi/2.)*np.arange(-1, 1+dtheta, dtheta)
    
    field_lines = []
    r = 1.
    for angle in angles:
        field_lines.append([r*np.sin(angle), 0., r*np.cos(angle)])
        field_lines[-1] = cx.GSMtoMAG(field_lines[-1], time, 'car', 'sph')

    cp.plot(time, 'p', 'xz', field_lines=field_lines, dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])

if False:
    cp.plot(time, 'xy', None, 'p', dx=0.5, dy=0.5, png=False)
    cp.plot(time, 'yz', None, 'p', dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])
    
    cp.plot(time, 'xz', None, 'p', dx=0.5, dy=0.5,
            xlims=[-30, 30], ylims=[-20, 20], zlims=[0, 10], png=False)
    
    r = 1.01
    mlat = 57.50
    mlong = 176.
    cp.plot(time, [r, mlat, mlong], None, 'p', dx=0.5, dy=0.5,
            xlims=[-30, 30], ylims=[-20, 20], png=False)
    
