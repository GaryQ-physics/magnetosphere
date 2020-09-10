"""
pip install joblib

"""

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import numpy as np
import biot_savart as bs
import cxtransform as cx
from units_and_constants import phys
from probe import probe
import util
from make_grid import make_grid, make_axes

def integrate(run, time_fname, mlat, mlon, para=True,
        xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.125, 
        tonpfile=False, returnAll=False, debug=True):
    """

    Returns (Btot)
    -------
    (3,) numpy array
    
        using kameleon, through probe, returns the result of biot savart 
        integration on regular grid given by a mesh of 1d arrays X,Y,Z 
        (specified by input parameters) for the field at x0 corresponding to 
        mlat, mlon on earth.

    
    Parameters
    ----------
    time : tuple/list/array
        the time at which this data is considered, used for coordinate transformations and to autogenerate filename for probe if no filename provided
    mlat : float
        magentic latitude of x0 (radius is 1 R_e).
    mlon : float
        magentic latitude of x0 
    filename : string, optional
        The default is None.
        if not None, this is the filename (full including path) passed to probe.
        if None, the filename is autogenerated by util.time2filename
    para : boolean, OPTIONAL
        if true, then computation run in parallel on multiple processors. Output
        returned is identical.
        The default is True.
    xlims : list/tuple, optional
        list/tuple of length 2
        X ranges from X[0]=xlims[0] to X[-1]=xlims[1]
        The default is None.
    (likewise ylims , zlims)
    d : float, optional
        if not None, it overrides dx, dy, and dz to all be d
        The default is 0.125
    dx : float, optional
        stepsize between points in X 
        The default is None.
    (likewise dy , dz)
    N : integer, optional
        if not None, it overrides Nx, Ny, and Nz to all be N
        The default is None.
    Nx : integer, optional
        the number of points in X 
        The default is None.
    (likewise Ny , Nz)
    L : float, optional
        if not None, it overrides xlims,ylims,and zlims to all be (-L,L)
        The default is None.
    fullVolume : boolean, OPTIONAL
        if True, it overides xlims, ylims and zlims to be the limits of the 
        full kameleon grid. 
        The default is False.
    fineVolume : bolean, OPTIONAL
        if True, it overides everything so that the grid is like the fine mesh 
        grid in the SWMF files except extend to be regular over a whole cube. 
        The default is False.
    print_output : boolean, optional
        if True: prints output on command line, and also if what is being run is 
        the standard test event (2003, 11, 20, 7, 0, 57.50, 176.00) then it appends
        results to a text file.
        The default is False.
    tolerance : floar, optional
        the tolerance to which which N and d type choices must be consistent.
        slight errors can be due to floating point arithmetic. 
        The default is 1e-13.
        Note if fineVolume=True, this is overidden to 0.

    """

    ax_list = make_axes(xlims, ylims, zlims, d)
    Nx = ax_list[0].size
    Ny = ax_list[1].size
    Nz = ax_list[2].size

    if Nx*Ny*Nz > 257**3:
        raise ValueError("number of points exceeds 257**3, potentially may run out of memory")

    if type(time_fname) == str:
        filename = os.path.split(time_fname)[1]
        time = util.CDFfilename2time(run, time_fname)
    else:
        time = time_fname
        filepath = util.time2CDFfilename(run, time_fname)

    util.dlfile(filepath, debug=True)

    x0 = cx.MAGtoGSM([1., mlat, mlon], time, 'sph', 'car')

    import tempfile
    os.system('rm ' + tempfile.gettempdir() + '/*dB_array_slice*')

    if para:
        #Gy, Gz = np.meshgrid(Y,Z)
        #Gy = Gy.flatten(order='C')
        #Gz = Gz.flatten(order='C')
        G_s = make_grid(ax_list, slices=True)

        def dBslice(i, debug=False):
            #Grid = np.column_stack([X[i]*np.ones(Gy.shape), Gy, Gz])
            J_kameleon = probe(filepath, G_s[i], var = ['jx','jy','jz'], library='kameleon')
            J = J_kameleon*(phys['muA']/phys['m']**2)
            if debug:
                print(G_s[i].shape)
                print(J.shape)

            dB_slice = bs.deltaB('dB', x0, G_s[i], J, V_char = dx*dy*dz)
            #deltaB_slice = np.sum(dB, axis=0)

            if tonpfile:
                npfname = tempfile.gettempdir() + '/dB_array_slice%d'%(i) + '.bin'
                #if os.path.exists(npfname):
                #    os.remove(npfname)
                dB_slice.tofile(npfname)

            return dB_slice

        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is not None and num_cores > X.size:
            num_cores = X.size
        print('Parallel processing {0:d} slices(s) using {1:d} cores'\
              .format(X.size, num_cores))
        dB_slices = Parallel(n_jobs=num_cores)(delayed(dBslice)(j) for j in range(X.size))

        G = np.column_stack(G_s).reshape((Nx*Ny*Nz,3)) #this is equivalent to having used slices=False
                                                        #see make_grid documentation, this 
        dB = np.column_stack(dB_slices).reshape((Nx*Ny*Nz,3)) # stitch together dB's in matching way

    else:
        G = make_grid(ax_list, slices=False)
        if debug:
            print('done G')

        J_kameleon = probe(filepath, G, var = ['jx','jy','jz'], library='kameleon')
        J = J_kameleon*(phys['muA']/phys['m']**2)
        if debug:
            print('done J')

        dB = bs.deltaB('dB', x0, G, J, V_char = d**3)
        if debug:
            print('done dB')

    if returnAll:
        return [dB, G, (Nx,Ny,Nz)]
    else:
        Btot = np.sum(dB, axis=0)
        return Btot

'''

    if print_output:
        print('Nx, Ny, Nz = {0:d}, {1:d}, {2:d}'.format(Nx,Ny,Nz))
        print('fullVolume = ' + str(fullVolume))
        print('fineVolume = ' + str(fineVolume))
        print('slices = ' + str(X.size))
        print('points in slice = ' + str(Y.size*Z.size))
        print('total points = ' + str(X.size*Y.size*Z.size))
        print('X[0], X[-1], dx = {0:f}, {1:f}, {2:f}'.format(X[0], X[-1], dx))
        print('Y[0], Y[-1], dy = {0:f}, {1:f}, {2:f}'.format(Y[0], Y[-1], dy))
        print('Z[0], Z[-1], dz = {0:f}, {1:f}, {2:f}'.format(Z[0], Z[-1], dz))
        print('para = ' + str(para))


    if print_output:
        print('Btot = \n' + str(Btot))
        print('Btot_norm = ' + str(np.linalg.norm(Btot)))

        string = '\n X[0], X[-1], dx; Y[0], Y[-1], dy; Z[0], Z[-1], dz; Btot[0], Btot[1], Btot[2], np.linalg.norm(Btot) = {0:f}, {1:f}, {2:f}, {3:f}, {4:f}, {5:f}, {6:f}, {7:f}, {8:f}, {9:f}, {10:f}, {11:f}, {12:f}'.format(X[0], X[-1], dx, Y[0], Y[-1], dy, Z[0], Z[-1], dz, Btot[0], Btot[1], Btot[2], np.linalg.norm(Btot))
        print(string)
        if tpad(time) == (2003, 11, 20, 7, 0, 0, 0) and (mlat, mlon) == (57.50, 176.00):
            datafname = conf['run_path_derived'] + 'biot_savart_kameleon_data_2003:11:20T07:00:00.txt'
            f = open(datafname,'a') # append only mode
            f.write(string)
            f.close()

    return Btot
'''

def dB_dV_slice(run, time_fname, mlat, mlon, u, v, U, returnAll=False):
    U1 = U[0]
    U2 = U[1]
    U3 = U[2]
    if isinstance(u, float):
        G = u*U1 + v*U2
    else:
        G = np.einsum('i,j->ij', u, U1) + np.einsum('i,j->ij', v, U2)
    

    if type(time_fname) == str:
        filename = os.path.split(time_fname)[1]
        time = util.CDFfilename2time(run, time_fname)
    else:
        time = time_fname
        filepath = util.time2CDFfilename(run, time_fname)

    util.dlfile(filepath, debug=True)

    x0 = cx.MAGtoGSM([1., mlat, mlon], time, 'sph', 'car')

    J_kameleon = probe(filepath, G, var = ['jx','jy','jz'], library='kameleon')
    J = J_kameleon*(phys['muA']/phys['m']**2)

    dB = bs.deltaB('dB', x0, G, J, V_char = 1.)

    if returnAll:
        return [dB, G, None]
    else:
        Btot = np.sum(dB, axis=0)
        return Btot
  


def toMAGLocalComponents(time, mlat, mlon, dB):
    if dB.shape == (3,):
        return toMAGLocalComponents(time, mlat, mlon, np.array([dB]))[0,:]

    time = np.array(time, dtype=int)

    if len(time.shape) == 1 and len(dB.shape) > 1:
        time = np.repeat([time], dB.shape[0], axis = 0)

    # time is Nx6 and dB is Nx3  mlat, mlon are numbers
    N = time.shape[0]
    print('time.shape == ' + str(time.shape))
    assert(time.shape[0] == N)
    assert(dB.shape == (N,3))
    station_pos = cx.MAGtoGSM(np.array([1., mlat, mlon]), time, 'sph', 'car')
    assert( station_pos.shape == (N, 3) ) # check

    Pole = cx.MAGtoGSM(np.array([0., 0., 1.]), time, 'car', 'car')

    U3 = station_pos
    U3norm = np.sqrt(U3[:,0]**2 + U3[:,1]**2 + U3[:,2]**2) #( not really needed since norm is 1 in these units where R_e=1, but just to be consistent in all units)
    divU3norm = 1./U3norm
    U3 = U3*divU3norm[:,np.newaxis]
    U1 = np.cross(Pole, U3)
    U1norm = np.sqrt(U1[:,0]**2 + U1[:,1]**2 + U1[:,2]**2)
    divU1norm = 1./U1norm
    #https://stackoverflow.com/questions/5795700/multiply-numpy-array-of-scalars-by-array-of-vectors
    U1 = U1*divU1norm[:,np.newaxis]
    U2 = np.cross(U3, U1)
    assert(U1.shape == (N, 3)) #check

    #print(U1)
    #print(U2)
    #print(U3)

    R = np.empty((N, 3, 3))
    R[:, :, 0] = U2
    R[:, :, 1] = U1
    R[:, :, 2] = -U3
    
    R = np.linalg.inv(R)
    
    #for i in range(N):
        #R[i,:,0] = U2[i,:]
        #R[i,:,1] = U1[i,:]
        #R[i,:,2] = -U3[i,:]
    #    M = np.column_stack([U2[i,:], U1[i,:], -U3[i,:]])
    #    R[i, :, :] = np.linalg.inv(M)

    dB_rot = np.einsum('ijk,ik->ij', R, dB)
    return dB_rot

