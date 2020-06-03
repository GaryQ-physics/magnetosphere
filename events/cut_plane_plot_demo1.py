import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import cut_plane_plot as cp

time = [2003, 11, 20, 7, 0, 0, 0]

cp.plot(time, 'xy', None, 'p', dx=0.5, dy=0.5)
cp.plot(time, 'xz', None, 'p', dx=0.5, dy=0.5, xlim=[-30, 30], ylim=[-20, 20])
cp.plot(time, 'yz', None, 'p', dx=0.5, dy=0.5, xlim=[-30, 30], ylim=[-20, 20])

cp.plot(time, 'xz', None, 'p', dx=0.5, dy=0.5,
        xlim=[-30, 30], ylim=[-20, 20], zlim=[0, 10])

r = 1.01
mlat = 57.50
mlong = 176.
cp.plot(time, [r, mlat, mlong], None, 'p', dx=0.5, dy=0.5,
        xlim=[-30, 30], ylim=[-20, 20], png=False)

