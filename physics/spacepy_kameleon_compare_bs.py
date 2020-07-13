'''
typical output:

opening: /home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070000-000.out
[ 0.16935409  3.49882983 -2.03944467]
4.053372119835279
[ 0.16935409  3.49882983 -2.03944467]
4.053372119895642
[ 0.47685288  3.48787187 -2.54728568]
4.345262153715399
[ 0.50112569  3.46145452 -2.60177497]
4.359131487480463

'''


import sys
import os
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import numpy as np

from units_and_constants import phys
import cxtransform as cx
import biot_savart as bs
from probe import probe
from util import time2filename

import spacepy.pybats.bats as bats


test_Event = np.array([2003, 11, 20, 7, 0, 57.50, 176.00])

time = test_Event[0:5]
mlat = test_Event[5]
mlon = test_Event[6]

filename = time2filename(time)[0:-4]
print('opening: ' + filename)
# read in the 3d magnetosphere
data3d = bats.Bats2d(filename)

# get the cell coordinates
x = data3d['x']
y = data3d['y']
z = data3d['z']
jx = data3d['jx']
jy = data3d['jy']
jz = data3d['jz']

x = np.array(x, dtype=float)
y = np.array(y, dtype=float)
z = np.array(z, dtype=float)
jx = np.array(jx, dtype=float)
jy = np.array(jy, dtype=float)
jz = np.array(jz, dtype=float)

D = 3.96875
xax = np.linspace(-D, D, 128)
yax = np.linspace(-D, D, 128)
zax = np.linspace(-D, D, 128)
dx = 0.0625
dV = dx**3

assert((D+D)/(128-1) == dx)
assert(xax[1] - xax[0] == dx)

Tr = np.all([-D<=x, x<=D, -D<=y, y<=D, -D<=z, z<=D], axis=0)
x_ = x[Tr]
y_ = y[Tr]
z_ = z[Tr]
jx_ = jx[Tr]
jy_ = jy[Tr]
jz_ = jz[Tr]

x0 = cx.MAGtoGSM([1., mlat, mlon], time, 'sph', 'car')

X = np.column_stack([x_, y_, z_])

J_sp = np.column_stack([jx_, jy_, jz_])
out_spacepy = bs.deltaB('deltaB', x0, X, J_sp*(phys['muA']/phys['m']**2), V_char = dV)
B_spacepy = np.linalg.norm(out_spacepy)

J_kam = probe(time, X, var=['jx', 'jy', 'jz'], usekV=True)
out_kam1 = bs.deltaB('deltaB', x0, X, J_kam*(phys['muA']/phys['m']**2), V_char = dV)

Grid = bs.make_grid(xax, yax, zax, 0, 0, 0)[0]
J_kam_grid = probe(time, Grid, var=['jx', 'jy', 'jz'], usekV=True)
out_kam2 = bs.deltaB('deltaB', x0, Grid, J_kam_grid*(phys['muA']/phys['m']**2), V_char = dV)

L = 4.
custGrid = bs.make_grid([-L, L], [-L, L], [-L, L], 0.1, 0.1, 0.1)[0]
J_kam_cust = probe(time, custGrid, var=['jx', 'jy', 'jz'], usekV=True)
out_kam3 = bs.deltaB('deltaB', x0, custGrid, J_kam_cust*(phys['muA']/phys['m']**2), V_char = 0.1**3)

print(out_spacepy)
print(np.linalg.norm(out_spacepy))
print(out_kam1)
print(np.linalg.norm(out_kam1))
print(out_kam2)
print(np.linalg.norm(out_kam2))
print(out_kam3)
print(np.linalg.norm(out_kam3))
