'''
typical output:

L = 4
    opening: /home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070000-000.out
    [ 0.16935409  3.49882983 -2.03944467]
    4.053372119835279
    [ 0.16935409  3.49882983 -2.03944467]
    4.053372119895642
    [ 0.47685288  3.48787187 -2.54728568]
    4.345262153715399
    [ 0.50112569  3.46145452 -2.60177497]
    4.359131487480463

    after importing biot_savart_kameleon:

    [ 0.16935409  3.49882983 -2.03944467]
    4.053372119835279
    [ 0.16935409  3.49882983 -2.03944467]
    4.053372119895642
    [ 0.47685288  3.48787187 -2.54728568]
    4.345262153715209
    [ 0.50112569  3.46145452 -2.60177497]
    4.359131487480463

L=8  ,  quick=True:
    [ 0.16935409  3.49882983 -2.03944467]
    4.053372119835279
    [ 0.16935409  3.49882983 -2.03944467]
    4.053372119895642
    [ 0.47685288  3.48787187 -2.54728568]
    4.345262153715209
    [-4.06456466 -3.77374938 10.48775325]
    11.864014436500856

L=8  ,  quick=False:
    [ 0.12519011  3.50014164 -1.96692771]
    4.016897890817517
    [ 0.12519011  3.50014164 -1.96692771]
    4.01689789087069
    [ 0.47685288  3.48787187 -2.54728568]
    4.345262153715209
    [-4.06456466 -3.77374938 10.48775325]
    11.864014436500856

L=4.  ,  quick=False:
    [ 0.12519011  3.50014164 -1.96692771]
    4.016897890817517
    [ 0.12519011  3.50014164 -1.96692771]
    4.01689789087069
    [ 0.47685288  3.48787187 -2.54728568]
    4.345262153715209
    [ 0.50112569  3.46145452 -2.60177497]
    4.359131487480397

NOTE: for this test_Event np.array([2003, 11, 20, 7, 0, 57.50, 176.00])
        the SWMF mag file has dBmhd = 4.581437951020489


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
import biot_savart_kameleon_interpolated_grid as bsk
import biot_savart_SWMF_native_finegrid as bss

import spacepy.pybats.bats as bats


test_Event = np.array([2003, 11, 20, 7, 0, 57.50, 176.00])

time = test_Event[0:5]
mlat = test_Event[5]
mlon = test_Event[6]

filename = time2filename(time).replace('.out.cdf','.out')
print('opening: ' + filename)
# read in the 3d magnetosphere
data3d = bats.Bats2d(filename)



#J_sp = np.column_stack([jx_, jy_, jz_])
#out_spacepy = bs.deltaB('deltaB', x0, X, J_sp*(phys['muA']/phys['m']**2), V_char = dV)
#B_spacepy = np.linalg.norm(out_spacepy)
out_spacepy = bss.run(data3d, time, mlat, mlon, quick=False)

dx = 0.0625
dV = dx**3
x0 = cx.MAGtoGSM([1., mlat, mlon], time, 'sph', 'car')
X = bss.run(data3d, time, mlat, mlon, returnX=True, quick=False)
J_kam = probe(time, X, var=['jx', 'jy', 'jz'], usekV=True)
out_kam1 = bs.deltaB('deltaB', x0, X, J_kam*(phys['muA']/phys['m']**2), V_char = dV)

#Grid = bs.make_grid(xax, yax, zax, 0, 0, 0)[0]
#J_kam_grid = probe(time, Grid, var=['jx', 'jy', 'jz'], usekV=True)
#out_kam2 = bs.deltaB('deltaB', x0, Grid, J_kam_grid*(phys['muA']/phys['m']**2), V_char = dV)
out_kam2 = bsk.integrate(time, mlat, mlon, para=True, spacepy_like=True, print_output=False)
# equivalent out_kam2 = bsk.integrate(time, mlat, mlon, para=True, L=3.96875, eps=0.0625, print_output=False, tolerance=0.)


L = 8.
#custGrid = bs.make_grid([-L, L], [-L, L], [-L, L], 0.1, 0.1, 0.1)[0]
#J_kam_cust = probe(time, custGrid, var=['jx', 'jy', 'jz'], usekV=True)
#out_kam3 = bs.deltaB('deltaB', x0, custGrid, J_kam_cust*(phys['muA']/phys['m']**2), V_char = 0.1**3)
out_kam3 = bsk.integrate(time, mlat, mlon, para=True, L=L, eps=0.1, print_output=False)

print(out_spacepy)
print(np.linalg.norm(out_spacepy))
print(out_kam1)
print(np.linalg.norm(out_kam1))
print(out_kam2)
print(np.linalg.norm(out_kam2))
print(out_kam3)
print(np.linalg.norm(out_kam3))
