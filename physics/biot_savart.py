import numpy as np
from units_and_constants import phys
from make_grid import make_grid, make_axes

import deltaB_subroutine


'''
import biot_savart as bs
import numpy as np
X0 = np.array([[0.,1,2],[2,2,2]])
X= 1.*np.arange(18).reshape(6,3)
J = 10*X
result = bs.deltaB('dB', X0, X, J)


result.shape

X0[0,:]
X0[1,:]


a=bs.deltaB('dB', X0[0,:], X, J)
b=bs.deltaB('dB', X0[1,:], X, J)
ap=bs.deltaB_old('dB', X0[0,:], X, J)
bp=bs.deltaB_old('dB', X0[1,:], X, J)
a==ap
b==bp



ap==result[0,:,:]
bp==result[1,:,:]
'''

def deltaB_wrap_Nx3(variable, X0, X, J, V_char = 1.):
    #X and J must be (N,3)
    #X0 must be (M,3)
    print('\n\n from deltaB_wrap_Nx3 \n\n')


def deltaB_wrap_3xN():
    pass


def deltaB(variable, X0, X, J, V_char = 1.):
    print('\n\n from deltaB \n\n')

    X=np.array(X)
    if X.shape == (3,):
        X=np.array([X])

    X0=np.array(X0)
    if X0.shape == (3,):
        X0=np.array([X0])

    X0 = np.repeat([X0], X.shape[0], axis=0)
    memloc = X0.__array_interface__['data'][0]
    X0=np.swapaxes(X0,0,1)
    print(memloc == X0.__array_interface__['data'][0])

    num_eval_pts = X0.shape[0]

    R = np.repeat([X], num_eval_pts, axis=0)
    del X
    R *= -1
    R += X0
    del X0

    dB = np.repeat([J], num_eval_pts, axis=0) #extended J right now, but will become dB after modifications

    rcut = 1.733*np.cbrt(V_char) # np.sqrt(3) == 1.7320508075688772
    Rcubed = np.einsum('ijk,ijk->ij',R,R)**1.5
    #divRcubed = 1./Rcubed
    #https://stackoverflow.com/questions/26248654/how-to-return-0-with-divide-by-zero/40022737
    divRcubed = np.divide(1., Rcubed, out=np.zeros_like(Rcubed), where=(Rcubed >= rcut**3))

    #dB = np.cross(dB, R)
    dB = dB.transpose()
    R = R.transpose()
    deltaB_subroutine.cross(dB, R)
    dB = dB.transpose()
    del R
    print('checkpoint')

    #dB = np.einsum('ijk,ij->ijk', dB, divRcubed)
    dB *= divRcubed[:,:,None]

    dB *= phys['mu0']/(4*np.pi)
    try:  #!!!!! better way?
        V_char.shape
        dB = np.einsum('i,ij->ij', V_char, dB)
    except:
        dB *= V_char

    if(variable=='dB'):
        if dB.shape[0]==1:
            dB = dB[0,:,:]
        return dB
    deltaB = np.sum(dB, axis=1)
    if(variable=='deltaB'):
        if deltaB.shape[0]==1:
            deltaB = deltaB[0,:]
        return deltaB
    return np.nan



def deltaB_old(variable, x0, X, J, V_char = 1.):
    print('\n\n from deltaB_old \n\n')
    X=np.array(X)
    if X.shape == (3,):
        X=np.array([X])

    x0=np.array(x0)
    

    X0 = np.repeat([x0], X.shape[0], axis=0)
    R = X0 - X

    rcut = 1.733*np.cbrt(V_char) # np.sqrt(3) == 1.7320508075688772
    Rcubed = (R[:,0]**2 + R[:,1]**2 + R[:,2]**2)**1.5
    #Rcubed = np.einstum('ij,ij->i',R,R)**1.5
    #divRcubed = 1./Rcubed
    #https://stackoverflow.com/questions/26248654/how-to-return-0-with-divide-by-zero/40022737
    divRcubed = np.divide(1., Rcubed, out=np.zeros_like(Rcubed), where=(Rcubed >= rcut**3))

    dB = (phys['mu0']/(4*np.pi))*( np.cross(J, R)*divRcubed[:,np.newaxis] ) #https://stackoverflow.com/questions/5795700/multiply-numpy-array-of-scalars-by-array-of-vectors
    #''''''''''''''''''''''''''''' np.einsum('ij,i->ij', np.cross(J, R), divRcubed)
    try:  #!!!!! better way?
        V_char.shape
        dB = np.einsum('i,ij->ij', V_char, dB)
    except:
        dB = V_char*dB

    if(variable=='dB'):
        return dB

    deltaB = np.sum(dB, axis=0)
    if(variable=='deltaB'):
        return deltaB

    return np.nan


def biot_savart_run(run, time, pts, regions, separate=False):
    import os
    import sys
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
    from config import conf
    import util
    from probe import probe

    filename = util.time2CDFfilename(run, time)

    if isinstance(regions, dict):
        regions = (regions,)

    ret = []
    for region in list(regions):
        ax_list = make_axes(region['xlims'], region['ylims'], region['zlims'], region['d'])
        G = make_grid(ax_list, slices=False)
        J = probe(filename, G, var=['jx','jy','jz'], library='kameleon')
        J *= phys['muA']/(phys['m']**2)
        result = deltaB('deltaB', pts, G, J, V_char = region['d']**3)
        ret.append(result.copy())
    ret = np.array(ret)
    if separate:
        return ret
    return np.sum(ret, axis=0)

if __name__ == '__main__':
    run = 'DIPTSUR2'
    time = (2019,9,2,6,30,0)

    #reg =  {'xlims': (-31.875, 46.875), 
    #        'ylims': (-31.875, 30.875),
    #        'zlims': (-31.875, 30.875), 
    #        'd': 0.25}

    #reg =  {'xlims': (-31.875, 31.875), 
    #        'ylims': (-31.875, 31.875),
    #        'zlims': (-31.875, 31.875), 
    #        'd': 0.25}

    reg = {'xlims': (-31.875, 31.875), 'd': 0.25, 'zlims': (-31.875, 30.875), 'ylims': (-31.875, 30.875)}

    ret = biot_savart_run(run, time, np.array([15., -1., -1.]), (reg,), separate=True)

    print(ret)
    print(ret.shape)

