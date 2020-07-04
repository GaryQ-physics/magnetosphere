# biot_savart_single

import numpy as np
from units_and_constants import phys


Rad = 0.5 # considered in units R_e
j_mag = 1. # considered in units muA/m^2
def J_kunits_single(X):
    if X[0]**2 + X[1]**2 < Rad**2:
        j = j_mag*np.array([0,0,1]) #in units of muA/m^2
    else:
        j = np.zeros((3,))
    return j

def deltaB_single(variable, X, X0, Npole, V_char = 1.):
    X=np.array(X)
    X0=np.array(X0)
    Npole=np.array(Npole)
    if np.dot(X0,X0)<1e-8: return np.nan
    J = J_kunits_single(X)*(phys['muA']/phys['m']**2)
    R = X0 - X # )*R_e

    #dB_dV = (mu0/(4*np.pi))*np.cross(J, R)/(np.linalg.norm(R)**3)
    #dBnT = V_char*dB_dV / nT
    dBnT = V_char*(phys['mu0']/(4*np.pi))*np.cross(J, R)/(np.linalg.norm(R)**3)
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

def B_EW_fromSingle(Xgrid, X0, Npole, dV_grid):
    Npole=np.array(Npole)

    a2 = np.cross(Npole, X0)
    a1 = np.cross(X0, a2)
    a1 = a1/np.linalg.norm(a1)
    a2 = a2/np.linalg.norm(a2)

    total = 0.
    for l in range(Xgrid.shape[0]):
        total = total + deltaB_single('dB_EW', Xgrid[l,:], X0, Npole, V_char = dV_grid)
    return total
