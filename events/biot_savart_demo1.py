# biot_savart_demo1

import numpy as np
from units_and_constants import phys
import biot_savart as bs
import biot_savart_single as bss

Npole = np.array([0., 0., 1.])

Rad = 0.5 # considered in units R_e
j_mag = 1. # considered in units muA/m^2
def J_kunits(X):
    X=np.array(X)
    if X.shape == (3,):
        X=np.array([X])

    Tr = X[:,0]**2 + X[:,1]**2 < Rad**2
    to_mult = Tr.astype(np.int) #https://www.python-course.eu/numpy_masking.php
    j = j_mag*np.array([0., 0., 1.]) # in units of muA/m^2
    j = np.repeat([j], X.shape[0], axis=0)
    j = j*to_mult[:,np.newaxis]
    return j


# analytic calculation:---------------
r_SI = 0.75*6.3781e+6

Rad_SI = Rad*6.3781e+6
J_SI = j_mag*1e-6
I_SI = J_SI*np.pi*Rad_SI**2 # 3.19501226e+7

Itot = np.pi*(Rad*phys['R_e'])**2*(j_mag*phys['muA']/phys['m']**2) # pi * (0.5*6.3781e+6 m)**2 * (1. muA/m^2) = 3.19501226e+13 muA = 3.19501226e+7 Amps
print('Itot= ', Itot)
print('I_SI= ', I_SI)
print(' ?=? 3.19501226')

B_SI = phys['mu0_SI']*I_SI/(2*np.pi*r_SI) # 1.33582614*1e-6
print('B_SI*1e+9= ', B_SI*1e+9) # 1335.82614

#----------------------

for i in range(4):
    mult = i+1.
    length = 20.*phys['R_e']

    dx = 0.05*phys['R_e']/mult
    dy = 0.05/mult
    dz = 0.1/mult
    xlims = [-1.1*Rad, 1.1*Rad]
    ylims = [-1.1*Rad, 1.1*Rad]
    zlims = [-0.5*length, +0.5*length]

    Xgrid = bs.make_grid(xlims, ylims, zlims, dx, dy, dz)
    X0 = np.array([0., 0.75, 0.])*phys['R_e']
    J = J_kunits(Xgrid)*(phys['muA']/phys['m']**2)

    print bs.B_EW(X0, Xgrid, J, Npole, dx*dy*dz)
    if i==0:
        print bss.B_EW_fromSingle(Xgrid, X0, Npole, dx*dy*dz)  # should be the same as line above but take longer
