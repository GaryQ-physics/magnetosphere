import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import cut_plane_plot as cp

'''
def plot(time, pos, plane_vs, parameter,
         xlims=[-10,10], ylims=[-10,10], zlims=None, 
         xticks=None, yticks=None, zticks=None,
         dx=0.1, dy=0.1,
         dpi=300, showplot=True,
         png=True, pngfile=None, debug=False):

    """
    plot(time, [r, mlat, mlong], None, 'p')
    plot(time, [GSMx,GSMy,GSMz], [v1, v2], 'p')
    """
'''

time = [2003, 11, 20, 7, 0, 0, 0]

cp.plot(time, 'xy', None, 'p', dx=0.5, dy=0.5, png=False)
cp.plot(time, 'xz', None, 'p', dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])
cp.plot(time, 'yz', None, 'p', dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])

cp.plot(time, 'xz', None, 'p', dx=0.5, dy=0.5,
        xlims=[-30, 30], ylims=[-20, 20], zlims=[0, 10], png=False)

r = 1.01
mlat = 57.50
mlong = 176.
cp.plot(time, [r, mlat, mlong], None, 'p', dx=0.5, dy=0.5,
        xlims=[-30, 30], ylims=[-20, 20], png=False)

