# biot_savart

import numpy as np
from units_and_constants import phys


def deltaB(variable, X0, X, J, V_char = 1.):
    X=np.array(X)
    if X.shape == (3,):
        X=np.array([X])

    X0=np.array(X0)
    

    X0 = np.repeat([X0], X.shape[0], axis=0)
    R = X0 - X

    Rcubed = (R[:,0]**2 + R[:,1]**2 + R[:,2]**2)**1.5
    divRcubed = 1./Rcubed

    dBnT = V_char*(phys['mu0']/(4*np.pi))*( np.cross(J, R)*divRcubed[:,np.newaxis] ) #https://stackoverflow.com/questions/5795700/multiply-numpy-array-of-scalars-by-array-of-vectors
    if(variable=='dB'):
        return dBnT

    deltaBnT = np.sum(dBnT, axis=0)
    if(variable=='deltaB'):
        return deltaBnT

    '''
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
    '''
    return np.nan

def B_EW(X0, X, J, Npole, dV_grid):
    Npole=np.array(Npole)

    a2 = np.cross(Npole, X0)
    a1 = np.cross(X0, a2)
    a1 = a1/np.linalg.norm(a1)
    a2 = a2/np.linalg.norm(a2)

    deltaBnT = deltaB('deltaB', X0, X, J, V_char=dV_grid)/phys['nT']
    return np.dot(deltaBnT,a2)

def make_grid(xlims, ylims, zlims, dx, dy, dz):
    no_origin = xlims[0] > 0. or xlims[1] < 0. or ylims[0] > 0. or ylims[1] < 0. or zlims[0] > 0. or zlims[1] < 0.
    if no_origin:
        print('WARNING: grid does not contain origin')
        X = np.arange(xlims[0], xlims[1]+dx, dx)
        Y = np.arange(ylims[0], ylims[1]+dy, dy)
        Z = np.arange(zlims[0], zlims[1]+dz, dz)
    else:
        X = np.concatenate([ -np.flip(np.delete(np.arange(0., -xlims[0]+dx, dx), 0)) , np.arange(0., xlims[1]+dx, dx) ])
        Y = np.concatenate([ -np.flip(np.delete(np.arange(0., -ylims[0]+dy, dy), 0)) , np.arange(0., ylims[1]+dy, dy) ])
        Z = np.concatenate([ -np.flip(np.delete(np.arange(0., -zlims[0]+dz, dz), 0)) , np.arange(0., zlims[1]+dz, dz) ])

    Nx= X.size
    Ny = Y.size
    Nz = Z.size
    B2, B3, B1 = np.meshgrid(Y, Z, X)
    #B1, B2, B3 = np.meshgrid(X, Y, Z) # seems more natural but doesnt work with vtk structured_grid format

    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')
    B3 = B3.flatten(order='C')
    Bgrid = np.column_stack((B1, B2, B3))

    return [Bgrid, Nx, Ny, Nz]
