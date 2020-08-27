import numpy as np
import time as tm

def make_grid(X, Y, Z, slices = False):

    if slices:
        Gy, Gz = np.meshgrid(Y,Z)
        Gy = Gy.flatten(order='C')
        Gz = Gz.flatten(order='C')

        ret = []
        for i in range(X.size):
            Grid = np.column_stack([X[i]*np.ones(Gy.shape), Gy, Gz])
            ret.append(Grid)
        return ret

    else: 
        Nx = X.size
        Ny = Y.size
        Nz = Z.size

        B2, B3, B1 = np.meshgrid(Y, Z, X)
        #B1, B2, B3 = np.meshgrid(X, Y, Z) # seems more natural but doesnt work with vtk structured_grid format
        B1 = B1.flatten(order='C')
        B2 = B2.flatten(order='C')
        B3 = B3.flatten(order='C')
        return np.column_stack((B1, B2, B3))



def paralelized_grids(X, Y, Z):

    Gy, Gz = np.meshgrid(Y,Z)
    Gy = Gy.flatten(order='C')
    Gz = Gz.flatten(order='C')

    ret = []
    for i in range(X.size):
        Grid = np.column_stack([X[i]*np.ones(Gy.shape), Gy, Gz])
        ret.append(Grid)

    return ret


def grid_in_biot_savart(X, Y, Z):
    Nx = X.size
    Ny = Y.size
    Nz = Z.size

    B2, B3, B1 = np.meshgrid(Y, Z, X)
    #B1, B2, B3 = np.meshgrid(X, Y, Z) # seems more natural but doesnt work with vtk structured_grid format
    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')
    B3 = B3.flatten(order='C')
    Bgrid = np.column_stack((B1, B2, B3))

    '''
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
    '''

    return Bgrid

n = 124 + 1
print(n)
#ar = np.arange(n)
ar = np.linspace(0,1,n)
X = ar
Y = ar
Z = ar

to = tm.time()
G1 = grid_in_biot_savart(X, Y, Z)
print(tm.time()-to)

to = tm.time()
G2l = paralelized_grids(X, Y, Z)
G2 = np.array(G2l)
G3 = np.column_stack(G2l)
G4 = G3.reshape((n**3,3))
print(tm.time()-to)

print(np.all(G4 == G1))

'''
print(G1.shape)
print(len(G2l))
print(G2l[0].shape)
print(G2.shape)
print(G3.shape)

print(G1)
print(G2l)
print(G2)
print(G3)
print(G3.flatten())
print(G1.flatten())
print(G3.reshape((n**3,3)))
print(G3.flatten() == G1.flatten())
print(G3.reshape((n**3,3)) == G1)
'''



