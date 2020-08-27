import numpy as np
import time as tm


def make_axes(xlims, ylims, zlims, d,
                dx=None, dy=None, dz=None,
                N=None,
                Nx=None, Ny=None, Nz=None,
                L=None,
                fullVolume=False, fineVolume=False, tolerance=1e-13):

    ###### make X, Y, and Z ###########################
    assert(not (fullVolume and spacepy))

    if fineVolume:
        L = 3.96875
        d = 0.0625
        N = 128
        assert((L+L)/(N-1) == d)
        assert((2*L)/(N-1) == d)
        tolerance = 0.
    if fullVolume:
        xlims = (-224., 32.)
        ylims = (-128., 128.)
        zlims = (-128., 128.)

    if L != None:
        xlims = (-L, L)
        ylims = (-L, L)
        zlims = (-L, L)

    if N != None:
        Nx = N
        Ny = N
        Nz = N

    if d != None:
        dx = d
        dy = d
        dz = d


    assert(xlims!=None and ylims!=None and zlims!=None)

    if Nx != None:
        X = np.linspace(xlims[0], xlims[1], Nx)
    elif dx != None:
        X = np.arange(xlims[0], xlims[1] + dx, dx)
        #X = np.arange(xlims[0], xlims[1], dx, endpoint=True) doesnt work
    else:
        assert(False)

    if Ny != None:
        Y = np.linspace(ylims[0], ylims[1], Ny)
    elif dy != None:
        Y = np.arange(ylims[0], ylims[1] + dy, dy)
    else:
        assert(False)

    if Nz != None:
        Z = np.linspace(zlims[0], zlims[1], Nz)
    elif dz != None:
        Z = np.arange(zlims[0], zlims[1] + dz, dz)
    else:
        assert(False)

    dx_check = X[1]-X[0]
    dy_check = Y[1]-Y[0]
    dz_check = Z[1]-Z[0]
    x_range_check = X[-1]-X[0]
    y_range_check = Y[-1]-Y[0]
    z_range_check = Z[-1]-Z[0]
    Nx_check = X.size
    Ny_check = Y.size
    Nz_check = Z.size

    assert(Nx_check == Nx or Nx == None)
    if dx != None:
        assert(np.abs(dx_check - dx) <= tolerance)
    assert(np.abs(xlims[1] - xlims[0] - x_range_check) <= tolerance)

    assert(Ny_check == Ny or Ny == None)
    if dy != None:
        assert(np.abs(dy_check - dy) <= tolerance or dy == None)
    assert(np.abs(ylims[1] - ylims[0] - y_range_check) <= tolerance)

    assert(Nz_check == Nz or Nz == None)
    if dz != None:
        assert(np.abs(dz_check - dz) <= tolerance or dz == None)
    assert(np.abs(zlims[1] - zlims[0] - z_range_check) <= tolerance)

    ax_list = [X, Y, Z]
    return ax_list


def make_grid(ax_list, slices=False):
    """
    note:   if G = make_grid(ax_list, slices=False)
            and G_s = make_grid(ax_list, slices=True)
            then G == np.column_stack(G_s).reshape((n**3,3))
    """
    X, Y, Z = ax_list

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

