import numpy as np


def dXds(X, s, sign, Bx_interp, By_interp, Bz_interp):

    if xlims[0]<X[0]<xlims[1] and ylims[0]<X[1]<ylims[1] and zlims[0]<X[2]<zlims[1]:
        #print('IN')
        B = np.array([Bx_interp(X)[0], By_interp(X)[0], Bz_interp(X)[0]])
        Bm = np.linalg.norm(B)
        if 1e-9 < Bm < 1e+7:
            return (sign/Bm)*B
    return [0., 0., 0.]


def fieldlines(run, time, mag, fieldvar='b', s_grid=None, max_iterations=100, debug=False):
    """Trace field lines from start points

    Parameters
    ----------
    time : list 
        or other type that is accepted by probe (cannot be 2d array)
    mag : numpy.array or list
        (Nx3) spherical mag coordinates of N points, which are the start points
        for N field lines
    s_grid : TYPE, optional
        DESCRIPTION. The default is np.arange(0., 200., 0.1)
    debug : TYPE, optional, boolean
        DESCRIPTION. The default is False.
    max_iterations : TYPE, optional, int
        stops integration after max_iterations amount, even if
        line goes out of bounds. The default is 100

    Returns
    -------
    ret : list of N arrays
        N field lines 

    """

    mag = np.array(mag)
    if len(mag.shape) == 1:
        mag = [mag]

    # Trace 3-D field line
    # TODO: Consider using
    # https://github.com/spacepy/spacepy/blob/master/spacepy/pybats/trace2d.py
    
    
    # make grid and interpolator on grid using probe ###################
    no_origin = xlims[0] > 0. or xlims[1] < 0. or ylims[0] > 0. or ylims[1] < 0. or zlims[0] > 0. or zlims[1] < 0.
    if no_origin:
        print('WARNING: grid does not contain origin')
        X = np.arange(xlims[0], xlims[1]+dx, dx)
        Y = np.arange(ylims[0], ylims[1]+dy, dy)
        Z = np.arange(zlims[0], zlims[1]+dz, dz)
    else:
        X = np.concatenate([ -np.flip(np.delete(np.arange(0., -xlims[0]+dx, dx), 0), 0) , np.arange(0., xlims[1]+dx, dx) ])
        Y = np.concatenate([ -np.flip(np.delete(np.arange(0., -ylims[0]+dy, dy), 0), 0) , np.arange(0., ylims[1]+dy, dy) ])
        Z = np.concatenate([ -np.flip(np.delete(np.arange(0., -zlims[0]+dz, dz), 0), 0) , np.arange(0., zlims[1]+dz, dz) ])
    Nx = X.size
    Ny = Y.size
    Nz = Z.size
    
    
    G2, G1, G3 = np.meshgrid(Y, X, Z) # different than in make_grid
    P = np.column_stack( (G1.flatten(), G2.flatten(), G3.flatten()) )
    filepath = time2CDFfilename(run, time, split=False)
    Fx = probe(filepath, P, var=fieldvar+'x').reshape(Nx, Ny, Nz)
    Fy = probe(filepath, P, var=fieldvar+'y').reshape(Nx, Ny, Nz)
    Fz = probe(filepath, P, var=fieldvar+'z').reshape(Nx, Ny, Nz)

    # https://stackoverflow.com/questions/21836067/interpolate-3d-volume-with-numpy-and-or-scipy
    Fx_interp = RegularGridInterpolator((X,Y,Z), Fx)
    Fy_interp = RegularGridInterpolator((X,Y,Z), Fy)
    Fz_interp = RegularGridInterpolator((X,Y,Z), Fz)
    ####################################################################
  
    
    if s_grid is None:
        s_grid = np.arange(0., 10., 0.1)
    # Trace field line for a total length of smax, and check if stop conditions
    # satified. If not satified, trace for another total length of smax.
    # Note that Python 3 version of integration library has stop function
    # that can be passed so this won't be needed.

    IC = cx.MAGtoGSM(mag, time[0:6], 'sph', 'car')
    ret = []
    linenum = 0
    for X0 in IC:
        if debug:
            print('linenum = ' + str(linenum))
        done = False
        solns = np.empty((0, 3)) # Combined solutions
        i = 0
        while not done:
            if debug:
                print('i = ' + str(i))
            soln = odeint(dXds, X0, s_grid, args=(-1, Fx_interp, Fy_interp, Fz_interp))
            R = soln[:, 0]**2+soln[:, 1]**2 + soln[:, 2]**2
            # define condition on the field line points
            # Find first location where soln steps out-of-bounds
            #tr = np.where( False == (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.) )        
            # Boolean array.


            tr = (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.)
            # RuntimeWarning: invalid value encountered in greater_equal


            # Indices where stop conditions satisfied
            tr_out = np.where(tr == False)
            if debug:
                print(tr)
            if tr_out[0].size > 0:
                # Stop condition found at least once. Use solution up to that point.s
                solns = np.vstack((solns, soln[0:tr_out[0][0] + 1, :]))
                done = True
            elif max_iterations == i + 1:
                solns = np.vstack((solns, soln))   # return soln   faster?
                done = True
            else:
                # New initial condition is stop point.
                X0 = soln[-1, :]
                # Append solution but exclude last value, which is the
                # new initial condition.
                solns = np.vstack((solns, soln[0:-1, :]))
            i = i + 1
        ret.append(solns)
        linenum += 1
            
    return ret
