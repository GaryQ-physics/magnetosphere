# biot_savart

import numpy as np
from units_and_constants import phys


def deltaB(variable, x0, X, J, V_char = 1.):
    X=np.array(X)
    if X.shape == (3,):
        X=np.array([X])

    x0=np.array(x0)
    

    X0 = np.repeat([x0], X.shape[0], axis=0)
    R = X0 - X

    Rcubed = (R[:,0]**2 + R[:,1]**2 + R[:,2]**2)**1.5
    divRcubed = 1./Rcubed

    dB = V_char*(phys['mu0']/(4*np.pi))*( np.cross(J, R)*divRcubed[:,np.newaxis] ) #https://stackoverflow.com/questions/5795700/multiply-numpy-array-of-scalars-by-array-of-vectors
    if(variable=='dB'):
        return dB

    deltaB = np.sum(dB, axis=0)
    if(variable=='deltaB'):
        return deltaB

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

def B_EW(x0, X, J, Npole, dV_grid):
    Npole=np.array(Npole)

    a2 = np.cross(Npole, x0)
    a1 = np.cross(x0, a2)
    a1 = a1/np.linalg.norm(a1)
    a2 = a2/np.linalg.norm(a2)
    print(list(a2))
    deltaBnT = deltaB('deltaB', x0, X, J, V_char=dV_grid)/phys['nT']
    return np.dot(deltaBnT,a2)

def make_grid(xlims, ylims, zlims, dx, dy, dz):
    if len(xlims) == 2:
        no_origin = xlims[0] > 0. or xlims[1] < 0. or ylims[0] > 0. or ylims[1] < 0. or zlims[0] > 0. or zlims[1] < 0.
        if no_origin:
            print('WARNING: grid does not contain origin')
            X = np.arange(xlims[0], xlims[1]+dx, dx)
            Y = np.arange(ylims[0], ylims[1]+dy, dy)
            Z = np.arange(zlims[0], zlims[1]+dz, dz)
        else:
            # need flip(a,0) for python 2.7.17 whereas flip(a,0) or flip(a) works in 2.7.18 and 3.7.4
            X = np.concatenate([ -np.flip(np.delete(np.arange(0., -xlims[0]+dx, dx), 0), 0) , np.arange(0., xlims[1]+dx, dx) ])
            Y = np.concatenate([ -np.flip(np.delete(np.arange(0., -ylims[0]+dy, dy), 0), 0) , np.arange(0., ylims[1]+dy, dy) ])
            Z = np.concatenate([ -np.flip(np.delete(np.arange(0., -zlims[0]+dz, dz), 0), 0) , np.arange(0., zlims[1]+dz, dz) ])
    else:
        X = np.array(xlims)
        Y = np.array(ylims)
        Z = np.array(zlims)
    Nx = X.size
    Ny = Y.size
    Nz = Z.size

    B2, B3, B1 = np.meshgrid(Y, Z, X)
    #B1, B2, B3 = np.meshgrid(X, Y, Z) # seems more natural but doesnt work with vtk structured_grid format
    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')
    B3 = B3.flatten(order='C')
    Bgrid = np.column_stack((B1, B2, B3))

    Xind = np.arange(0, Nx, 1)
    Yind = np.arange(0, Ny, 1)
    Zind = np.arange(0, Nz, 1)
    C2, C3, C1 = np.meshgrid(Yind, Zind, Xind)
    C1 = C1.flatten(order='C')
    C2 = C2.flatten(order='C')
    C3 = C3.flatten(order='C')
    inds = np.column_stack((C1, C2, C3))

    Xc = np.arange(0, Nx-1, 1)
    Yc = np.arange(0, Ny-1, 1)
    Zc = np.arange(0, Nz-1, 1)
    D2, D3, D1 = np.meshgrid(Yc, Zc, Xc)
    D1 = D1.flatten(order='C')
    D2 = D2.flatten(order='C')
    D3 = D3.flatten(order='C')
    cell_inds = np.column_stack((D1, D2, D3))

    return [Bgrid, Nx, Ny, Nz, cell_inds, inds]
