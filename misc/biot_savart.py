# biot_savart

import numpy as np

# units and constants---------------
hr = 1.
muA = 1.
R_e = 1.
nT = 1.

deg = (np.pi/180.)
amin = deg/60.
minn = hr/60.
s = minn/60.

A = 1e+6*muA
Tesla = 1e+9*nT
m = R_e/6.3781e+6  # R_e == 6.3781e+6*m  (note m < 10^-6 Re)
#kg = Tesla*A*s**2

mu0_SI = 1.25663706212e-6 # almost 4*pi*1e-7
#mu0 = mu0_SI*kg*m/((s**2)*(A**2))
mu0 = 1.970237314e-10*nT*R_e/muA
#--------------------

Npole = np.array([0., 0., 1.])
Rad = 0.5 # considered in units R_e
j_mag = 1. # considered in units muA/m^2

# analytic calculation:---------------
r_SI = 0.75*6.3781e+6

Rad_SI = 0.5*6.3781e+6
J_SI = j_mag*1e-6
I_SI = J_SI*np.pi*Rad_SI**2 # 3.19501226e+7

Itot = np.pi*(Rad*R_e)**2*(j_mag*muA/m**2) # pi * (0.5*6.3781e+6 m)**2 * (1. muA/m^2) = 3.19501226e+13 muA = 3.19501226e+7 Amps
print('Itot= ', Itot)
print('I_SI= ', I_SI)
print(' ?=? 3.19501226')

B_SI = mu0_SI*I_SI/(2*np.pi*r_SI) # 1.33582614*1e-6
print('B_SI*1e+9= ', B_SI*1e+9) # 1335.82614

#----------------------

def J_kunits(X):
    X=np.array(X)
    Tr = X[:,0]**2 + X[:,1]**2 < Rad**2
    to_mult = Tr.astype(np.int) #https://www.python-course.eu/numpy_masking.php
    j = j_mag*np.array([0., 0., 1.]) # in units of muA/m^2
    j = np.repeat([j], X.shape[0], axis=0)
    j = j*to_mult[:,np.newaxis]
    return j

def J_kunits_single(X):
    if X[0]**2 + X[1]**2 < Rad**2:
        j = j_mag*np.array([0,0,1]) #in units of muA/m^2
    else:
        j = np.zeros((3,))
    return j

def deltaB(variable, X, X0, Npole, V_char = 1.):
    X=np.array(X)
    if X.shape == (3,):
        X=np.array([X])

    X0=np.array(X0)
    Npole=np.array(Npole)

    a2 = np.cross(Npole, X0)
    a1 = np.cross(X0, a2)

    a1 = a1/np.linalg.norm(a1)
    a2 = a2/np.linalg.norm(a2)
    #a1 = np.repeat([a1], X.shape[0], axis=0)
    #a2 = np.repeat([a2], X.shape[0], axis=0)

    X0 = np.repeat([X0], X.shape[0], axis=0)
    J = J_kunits(X)*(muA/m**2)
    R = X0 - X # )*R_e

    Rcubed = (R[:,0]**2 + R[:,1]**2 + R[:,2]**2)**1.5
    divRcubed = 1./Rcubed

    dBnT = V_char*(mu0/(4*np.pi))*( np.cross(J, R)*divRcubed[:,np.newaxis] ) #https://stackoverflow.com/questions/5795700/multiply-numpy-array-of-scalars-by-array-of-vectors
    if(variable=='dB_seperately'):
        return dBnT
    deltaBnT = np.sum(dBnT, axis=0)
    if(variable=='deltaB'):
        return deltaBnT
    if(variable=='deltaB_mag'):
        return np.sqrt(deltaBnT[:,0]**2 + deltaBnT[:,1]**2 + deltaBnT[:,2]**2)
    if(variable=='deltaB_EW'):
        # east west direction (east positive)
        #return np.einsum('ij,ij->i', deltaBnT, a2) #https://stackoverflow.com/questions/15616742/vectorized-way-of-calculating-row-wise-dot-product-two-matrices-with-scipy
        return np.dot(deltaBnT,a2)
    if(variable=='deltaB_NS'):
        # north south direction (north positive)
        return np.einsum('ij,ij->i', deltaBnT, a1)
    if(variable=='deltaBx'):
        return deltaBnT[:,0]
    return np.nan

def deltaB_single(variable, X, X0, Npole, V_char = 1.):
    X=np.array(X)
    X0=np.array(X0)
    Npole=np.array(Npole)
    if np.dot(X0,X0)<1e-8: return np.nan
    J = J_kunits_single(X)*(muA/m**2)
    R = X0 - X # )*R_e

    #dB_dV = (mu0/(4*np.pi))*np.cross(J, R)/(np.linalg.norm(R)**3)
    #dBnT = V_char*dB_dV / nT
    dBnT = V_char*(mu0/(4*np.pi))*np.cross(J, R)/(np.linalg.norm(R)**3)
    if(variable=='dB'):
        return dBnT
    if(variable=='dB_mag'):
        return np.linalg.norm(dBnT)
    a2 = np.cross(Npole, X0)
    if(variable=='dB_EW'):
        return np.dot(dBnT, a2)/np.linalg.norm(a2) # east west direction (east positive)
    a1 = np.cross(X0, a2)
    if(variable=='dB_NS'):
        return np.dot(dBnT, a1)/np.linalg.norm(a1) # north south direction (north positive)
    if(variable=='dBx'):
        return dBnT[0]


if False: # test that the deltaB_single function and variety of inputs on the deltaB function agree 
    print deltaB('dB_seperately', [[0,0.5,0], [0,.499,0], [0,0,1]], [0,2,0], Npole)

    print deltaB('dB_seperately', [[0,0.5,0]], [0,2,0], Npole)
    print deltaB('deltaB', [[0,0.5,0]], [0,2,0], Npole)
    print deltaB('dB_seperately', [0,0.5,0], [0,2,0], Npole)
    print deltaB('deltaB', [0,0.5,0], [0,2,0], Npole)

    print deltaB('dB_seperately', [[0,.499,0]], [0,2,0], Npole)
    print deltaB('deltaB', [[0,.499,0]], [0,2,0], Npole)
    print deltaB('dB_seperately', [0,.499,0], [0,2,0], Npole)
    print deltaB('deltaB', [0,.499,0], [0,2,0], Npole)

    print deltaB('dB_seperately', [[0,0,1]], [0,2,0], Npole)
    print deltaB('deltaB', [[0,0,1]], [0,2,0], Npole)

    print deltaB_single('dB', [0,0.5,0], [0,2,0], Npole)
    print deltaB_single('dB', [0,.499,0], [0,2,0], Npole)
    print deltaB_single('dB', [0,0,1], [0,2,0], Npole)

    print deltaB_single('dB', [0,0.5,0], [0,2,0], Npole) + deltaB_single('dB', [0,.499,0], [0,2,0], Npole) + deltaB_single('dB', [0,0,1], [0,2,0], Npole)
    print deltaB('deltaB', [[0,0.5,0], [0,.499,0], [0,0,1]], [0,2,0], Npole)


def B_EW(X0, Npole, mult=1, length=10.):
    dx = 0.05/mult
    dy = 0.05/mult
    dz = 0.1/mult
    X = np.arange(-Rad - 0.1, Rad + 0.1 + dx, dx)
    Nx = X.size # 24
    Y = np.arange(-Rad - 0.1, Rad + 0.1 + dy, dy)
    Ny = Y.size
    Z = np.arange(-0.5*length, 0.5*length + dz, dz)
    Nz = Z.size # 100
    B2, B3, B1 = np.meshgrid(Y, Z, X)
    #B1, B2, B3 = np.meshgrid(X, Y, Z)

    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')
    B3 = B3.flatten(order='C')
    Bgrid = np.column_stack((B1, B2, B3))

    total = deltaB('deltaB_EW', Bgrid, X0, Npole, V_char = dx*dy*dz) # dx*dy*dz*R_e**3
    return total

def B_EW_fromSingle(X0, Npole, mult=1, length=10.):
    dx = 0.05/mult
    dy = 0.05/mult
    dz = 0.1/mult
    X = np.arange(-Rad - 0.1, Rad + 0.1 + dx, dx)
    Nx = X.size # 24
    Y = np.arange(-Rad - 0.1, Rad + 0.1 + dy, dy)
    Ny = Y.size
    Z = np.arange(-0.5*length, 0.5*length + dz, dz)
    Nz = Z.size # 100
    B2, B3, B1 = np.meshgrid(Y, Z, X)
    #B1, B2, B3 = np.meshgrid(X, Y, Z)

    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')
    B3 = B3.flatten(order='C')
    Bgrid = np.column_stack((B1, B2, B3))

    total = 0.
    for l in range(B1.size):
        total = total + deltaB_single('dB_EW', Bgrid[l,:], X0, Npole, V_char = dx*dy*dz) # dx*dy*dz*R_e**3
    return total

if False:
    print B_EW([0., 0.75, 0.], Npole, mult=i+1)
    print B_EW_fromSingle([0., 0.75, 0.], Npole)  # should be the same as line above but take longer

for i in range(8):
    print B_EW([0., 0.75, 0.], Npole, mult=i+1)
print ('---')
for i in range(4):
    print B_EW([0., 0.75, 0.], Npole, mult=4, length=10.*(i+1))
