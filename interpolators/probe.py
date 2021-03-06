import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

#from util import time2filename, filemeta
import util
from units_and_constants import phys

'''
import probe as p
import numpy as np
di = np.linspace(0,15,151)
X = np.zeros((151,3))
X[:,0] = di
X[:,2] = di
p.J_analytic(X)
p.B_analytic(X)
'''

def J_analytic(X):
    Mhat = np.array([0,0,1.])

    if len(X.shape)==1:
        X = np.array([X])

    '''
    Jchar = 10. # thought of in kameleon units (muA/m^2)
    a = 5. #X and a in same units (R_e, a base unit)
    def jrad(r):
        return Jchar*np.exp(-a*r)
    '''

    rCurrents = 1.5
    Jchar = 1. # thought of in kameleon units (muA/m^2)
    a = 3. #X and a in same units (R_e, a base unit)
    def jrad(r):
        divr = np.divide(1., r, out=np.zeros_like(r), where=(r >= rCurrents))
        return Jchar*a**5*divr**5

    R = np.sqrt(np.einsum('ij,ij->i', X, X))
    divR = np.divide(1., R, out=np.zeros_like(R), where=(R != 0.))
    J = np.einsum('i,ij->ij', (jrad(R) * divR), np.cross(Mhat, X) )

    if J.shape[0]==1:
        J = J[0,:]
    return J #return thought in kameleon units since Jchar is

def B_analytic(X):
    Mhat = np.array([0,0,1.])

    if len(X.shape)==1:
        X = np.array([X])

    '''
    Jchar = 10.
    a = 5.
    Jchar_K = Jchar * ( phys['muA']/(phys['m']**2) )
    def Bc(r):
        #integrate mu0/3 *( Jrad(x) ) from r to infinity
        return (phys['mu0']/3.)*Jchar_K*(1./a)*np.exp(-a*r)
    def Md(r):
        #integrate mu0/3 *( x^3 * Jrad(x) ) from 0 to r
        return (phys['mu0']/3.)*Jchar_K*( 6/a**4 - (a*r*(a*r*(a*r+3)+6)+6)*np.exp(-a*r)/a**4 )
    '''

    rCurrents = 1.5
    Jchar = 1.
    a = 3.
    Jchar_K = Jchar * ( phys['muA']/(phys['m']**2) )
    def Bc(r):
        #integrate mu0/3 *( Jrad(x) ) from r to infinity
        divr = np.divide(1., r, out=np.zeros_like(r), where=(r >= rCurrents))
        ret = (phys['mu0']/3.)*Jchar_K*(a**5/4) * divr**4
        ret[divr==0.] = (phys['mu0']/3.)*Jchar_K*(a**5/4) / (rCurrents**4)
        return ret
    def Md(r):
        #integrate mu0/3 *( x^3 * Jrad(x) ) from 0 to r (or rather rCurrents to r since take J=0 below that
        divr = np.divide(1., r, out=np.zeros_like(r), where=(r >= rCurrents))
        ret = (phys['mu0']/3.)*Jchar_K*a**5 *(1./rCurrents - divr)
        ret[divr==0.] = 0.
        return ret

    # B is sum of two terms, on from integrating dipole field from shells inside (using Md),
    # and other from integrating constant field from shells outside (using Bc)
    R = np.sqrt(np.einsum('ij,ij->i', X, X))
    divR = np.divide(1., R, out=np.zeros_like(R), where=(R > 0.))

    #print(X)
    #print(R)
    #rint(divR)
    #print(Md(R))
    #print(Bc(R))

    B = (3* (Md(R)*np.einsum('j,ij',Mhat, X)*divR**5)[:,None] * X - Mhat*(Md(R)*divR**3)[:,None] ) \
      + (2*Mhat*Bc(R)[:,None])

    if B.shape[0]==1:
        B = B[0,:]
    return B

'''
import numpy as np
import probe as p
X = 10.*np.random.rand(18).reshape((6,3))
p.B_analytic_loop(X[0,:])
p.B_analytic(X[0,:])
p.B_analytic_loop(X)
p.B_analytic(X)
np.allclose(p.B_analytic_loop(X), p.B_analytic(X))
np.sqrt(np.einsum('ij,ij->i', X, X))
'''

def B_analytic_loop(X):
    if X.shape != (3,):
        toret = np.nan*np.empty(X.shape)
        for i in range(X.shape[0]):
            toret[i,:] = B_analytic_loop(X[i,:])
        return toret


    Mhat = np.array([0,0,1.])
    rCurrents = 1.5
    Jchar = 1.
    a = 3.
    Jchar_K = Jchar * ( phys['muA']/(phys['m']**2) )
    def Bc(r):
        #integrate mu0/3 *( Jrad(x) ) from r to infinity
        if r >=rCurrents:
            return (phys['mu0']/3.)*Jchar_K*((a**5)/4.) / (r**4)
        else:
            return (phys['mu0']/3.)*Jchar_K*((a**5)/4.) / (rCurrents**4)

    def Md(r):
        #integrate mu0/3 *( x^3 * Jrad(x) ) from 0 to r (or rather rCurrents to r since take J=0 below that
        if r >=rCurrents:
            return (phys['mu0']/3.)*Jchar_K*a**5 *(1./rCurrents - 1/r)
        else:
            return 0.

    R = np.linalg.norm(X)

    B = ( 3*(Md(R)* np.dot(Mhat,X) / (R**5)) * X - (Md(R)/ (R**3)) * Mhat ) \
      + ( 2*Bc(R)*Mhat )

    return B


def J_analytic2(X):
    """
    X
        Nx3 array of N 3 vectors, assumed in units of R_e
    """

    M = np.array([0,0,1000.]) * phys['nT']*phys['R_e']**3

    amin = 2. * phys['R_e']
    amax = 3. * phys['R_e']
    psi = M * 15./( phys['mu0']*(amax**5-amin**5) ) # psi = rho*omega (charge density * angular velocity)
    if False:
        B_inside = M*5*(amax**2-amin**2)/(amax**5-amin**5) # == (mu0/3)*(amax**2-amin**2)*psi
        print(B_inside) # 1000.*5*(3**2-2**2)/(3**5-2**5) --> 118.48341232227489

    Rsq = np.einsum('ij,ij->i', X, X)
    Tr = np.logical_and(amin**2 < Rsq, Rsq < amax**2)
    to_mult = Tr.astype(np.int) #https://www.python-course.eu/numpy_masking.php
    j = np.cross(X, psi) # as pure number this is in units of muA/(R_e**2) since those are base units of phys
    j = np.einsum('i,ij->ij', to_mult, j) # 

    j_kameleonUnits = j/( phys['muA']/(phys['m']**2) )

    #print('HELLO\n\n\n\n\n\n')
    #print(np.max(j_kameleonUnits))
    #print(np.min(j_kameleonUnits))
    #print(j_kameleonUnits)

    return j_kameleonUnits


def probe(filename, P, var=None, debug=False, dictionary=False, library='kameleonV'):
    """
    library = 'kameleonV', 'kameleon', 'pycdf'
    """
    P = np.array(P)
    if P.shape == (3, ):
        P = np.array([P])

    assert(filename[0] == '/')
    TESTANALYTIC = 'TESTANALYTIC' in filename
    if debug: print('TESTANALYTIC = %s'%(str(TESTANALYTIC)))


    if not os.path.exists(filename):# and not TESTANALYTIC:
        raise ValueError('Not found: ' + filename)

    ################### import apropriate file for library
    if library == 'kameleonV':
        assert(P.shape[1] == 3)
        sys.path.append(conf['interpolator'] + 'kameleonV-0.2.3/')
        import kameleonV
    if library == 'kameleon':
        assert(P.shape[1] == 3)
        sys.path.append(conf['interpolator'] + 'kameleon/lib/python2.7/site-packages/ccmc/')
        if debug: print('before import')
        import _CCMC as ccmc
        if debug: print('after import')
    if library == 'pycdf':
        assert(P.shape[1] == 2)
        sys.path.append(conf['interpolator'] + 'pycdf_with_scipy')
        import ionosphere_interpolator as ii
    ###################


    ################### define interpolate(variable, Q) funtion within this probe
    if TESTANALYTIC:
        def interpolate(variable, Q):
            if variable == 'jx':
                return J_analytic(Q)[:, 0]
            if variable == 'jy':
                return J_analytic(Q)[:, 1]
            if variable == 'jz':
                return J_analytic(Q)[:, 2]
            if variable in ['bx', 'b1x']:
                return B_analytic(Q)[:, 0]
            if variable in ['by', 'b1y']:
                return B_analytic(Q)[:, 1]
            if variable in ['bz', 'b1z']:
                return B_analytic(Q)[:, 2]

    elif library == 'kameleonV':
        def interpolate(variable, Q):
            return kameleonV.interpolate(filename, Q[:,0], Q[:,1], Q[:,2], variable)

    elif library == 'kameleon':
        if debug: print('bef int')
        kameleon = ccmc.Kameleon()
        kameleon.open(filename)
        interpolator = kameleon.createNewInterpolator()

        def interpolate(variable, Q):
            kameleon.loadVariable(variable)
            arr = np.nan*np.empty((Q.shape[0],))
            for k in range(Q.shape[0]):
                arr[k] = interpolator.interpolate(variable, Q[k,0], Q[k,1], Q[k,2])
            return arr
        if debug: print('aft int')

    elif library == 'pycdf':
        def interpolate(variable, Q):
            return ii.interpolate(Q[:,0], Q[:,1], variable, filename)
    ###################

    if var is None:
        meta = util.filemeta(filename)
        # Get data for all parameters, store in dictionary
        ret = {}
        for key in meta['parameters']:
            ret[key] = interpolate(key, P)

    elif type(var) == str:
        ret = interpolate(var, P)
        if ret.size == 1:
            ret = ret[0]

    else:
        if type(var) != list and type(var) != tuple:
            raise ValueError('var must be None, str, list, or tuple')
        if dictionary:
            ret = {}
        else:
            ret = np.nan*np.empty((P.shape[0], len(var)))
        i = 0
        for v in var:
            if type(v) != str:
                raise ValueError('var must be a str, or a list/tuple of strs')              
            arr = interpolate(v, P)
            if arr.size == 1:
                arr = arr[0]
            if dictionary:
                ret[v] = arr
            else:
                ret[:,i] = arr
                i = i + 1
        if not dictionary:
            if P.shape[0] == 1:
                ret = ret.flatten()

    if library == 'kameleon' and not TESTANALYTIC:
        kameleon.close()
    if debug: print('DONE PROBING')
    return ret


def GetRunData(run, time, P, var):
    if var=='j':
        var = ['jx','jy','jz']
    elif var=='b':
        var = ['bx','by','bz']
    elif var=='b1':
        var = ['b1x','b1y','b1z']

    '''
    if 'TESTANALYTIC' == run:
        TESTANALYTIC = True
    else:
        TESTANALYTIC = False
    '''

    filename = util.time2CDFfilename(run, time)
    return probe(filename, P, var=var, library='kameleon')

