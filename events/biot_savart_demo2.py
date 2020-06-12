
# biot_savart_demo2

import numpy as np
from units_and_constants import phys
import biot_savart as bs
import biot_savart_single as bss

run_seperately = True

Npole = np.array([0., 0., 1.])
Vect = np.array([0., 0., 1.])

Rad = 0.5 # considered in units R_e
j_mag = 1. # considered in units muA/m^2
def J_kunits(X):
    X=np.array(X)
    if X.shape == (3,):
        X=np.array([X])

    Tr = np.einsum('ij,ij->i', X, X) < Rad**2
    to_mult = Tr.astype(np.int) #https://www.python-course.eu/numpy_masking.php
    j = np.cross(X, Vect) # in units of muA/m^2
    j = j*to_mult[:,np.newaxis]
    #print X[1300,:]
    #print Tr[1300]
    #print j[1300,:]
    return j


if run_seperately:
    # analytic calculation:---------------
    '''
    Vect_SI must be st:
        1e+6 * Vect_SI * (x/m)  =  Vect * (x/Re) =  J/(muA/m^2)
    '''
    Vect_SI = (1e-6*phys['m']/phys['R_e']) * Vect  # get from muA to A (already per m^2)
    Rad_SI = Rad * phys['R_e']/phys['m'] # Rad*6.3781e+6
    print(Rad_SI, Rad*6.3781e+6)
    moment = (4.*np.pi/15.)*(Rad_SI**5)*Vect_SI # dipole moment in SI units
    print(moment)

    R_SI = np.array([0., 0.75, 0.])*6.3781e+6
    r_SI = np.sqrt(np.dot(R_SI,R_SI)) # 0.75*6.3781e+6
    print(r_SI, 0.75*6.3781e+6)
    B_SI = (phys['mu0_SI']/(4.*np.pi))*( (3./r_SI**5)*np.dot(moment,R_SI)*R_SI - (1./r_SI**3)*moment )

    print('B_SI*1e+9= ', B_SI*1e+9) # B in nanotesla 
    #----------------------

    for i in range(13):
        mult = i+1.

        dx = 0.05*phys['R_e']/mult
        dy = 0.05/mult
        dz = 0.1/mult
        xlims = [-1.1*Rad, 1.1*Rad]
        ylims = [-1.1*Rad, 1.1*Rad]
        zlims = [-1.1*Rad, 1.1*Rad]

        Xgrid = bs.make_grid(xlims, ylims, zlims, dx, dy, dz)[0]
        X0 = np.array([0., 0.75, 0.])*phys['R_e']
        J = J_kunits(Xgrid)*(phys['muA']/phys['m']**2)

        #print J[1300,:]
        print bs.deltaB('deltaB', X0, Xgrid, J, V_char=dx*dy*dz)
        print bs.B_EW(X0, Xgrid, J, Npole, dx*dy*dz)
