# using Luhmann 1978 [https://doi.org/10.1029/JA084iA08p04405]
import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../jitFORTRAN')

from config import conf
import jitFORTRAN
from units_and_constants import phys

B0_FORTRAN = '''!B0_FORTRAN
subroutine B0_FORTRAN(M,X,RET,N)
    implicit none
    !will pass the transposed np array (which is in fortran order) having shape (3, N, M)

    integer, intent(in) :: N
    real(4), intent(in), dimension(0:2) :: M
!f2py intent(in) :: M
    real(4), intent(in), dimension(0:2, 0:N-1) :: X
!f2py intent(in) :: X
    real(4), intent(inout), dimension(0:2, 0:N-1) :: RET
!f2py intent(in,out) :: RET

    integer :: i
    real(4) :: divr2

    do i=0,N-1
        divr2 = 1./DOT_PRODUCT(X(:,i), X(:,i))
        RET(:,i) = ( 3.*DOT_PRODUCT(M, X(:,i))*divr2**2.5 )* X(:,i) - (divr2**1.5) * M
    end do

end subroutine B0_FORTRAN
'''
B0_F = jitFORTRAN.Fortran_Subroutine(B0_FORTRAN, 'B0_FORTRAN')

def B0(X): # dipole field
    M = np.array([0,0,3.12e+4], dtype=np.float32)#, dtype=np.float32) # "dipole moment"(not really) in  nT * R_e**3  # https://en.wikipedia.org/wiki/Dipole_model_of_the_Earth%27s_magnetic_field

    if X.shape == (3,):
        divr = 1./np.linalg.norm(X)
        return ( 3.*np.dot(M,X)*divr**5 )* X  -  (divr**3) * M

    if not X.dtype==np.float32:
        X = np.array(X, dtype=np.float32)

    ret = np.empty(X.shape, dtype=np.float32)
    print(np.isfortran(ret))
    N = X.shape[0]
    X = X.transpose()
    ret = ret.transpose()
    print(np.isfortran(ret))
    B0_F.execute(M, X, ret, N)
    X = X.transpose()
    ret = ret.transpose()
    print(np.isfortran(ret))
    return ret

def B1(X): # other part of field in L1978
    B_T = 100. # in nT
    delta = 3. # in R_e

    if X.shape == (3,):
        print('using one')
        delta = np.float32(delta)
        return B_T*np.tanh(X[2]/delta) * np.array([1.,0,0])

    if not X.dtype==np.float32:
        X = np.array(X, dtype=np.float32)

    ret = np.zeros(X.shape, dtype=np.float32)
    ret[:, 0] = B_T*np.tanh(X[:,2]/delta)
    return ret

def B(X):
    return B0(X) + B1(X)

def JA(X): # analytic J = curl B
    B_T = 100. # in nT
    delta = 3. # in R_e
    coef = ( B_T/(delta*phys['mu0']) )/( phys['muA']/(phys['m']**2) )

    # regardless of X.dtype or X.shape, output will be float32
    # however, for X.dtype=float32, if delta is float64 then 
    #               X[2]/delta will be float62 
    #       while X[:,2]/delta will be float32      (the term in cosh)

    if X.shape == (3,):
        print('using one')
        delta = np.float32(delta)
        return coef*np.cosh(X[2]/delta)**-2 * np.array([0,1.,0], dtype=np.float32)

    if not X.dtype==np.float32:
        X = np.array(X, dtype=np.float32)

    ret = np.zeros(X.shape, dtype=np.float32)
    ret[:, 1] = coef*np.cosh(X[:,2]/delta)**-2
    return ret


if __name__=='__main__':
    if False:
        P = 1.*np.random.rand(21).reshape((7,3))
        P = np.array(P, dtype=np.float32)
        print(P)
        vals = B0(P)
        print(vals)
        for i in range(7):
            print(JA(P[i,:]).dtype)
            print(vals[i,:] == B0(P[i,:]))
            print( ( vals[i,:] - B0(P[i,:]) )/vals[i,:] )

    if sys.argv[1] == 'pyth3':
        pyth3 = True
    elif sys.argv[1] == 'pyth2':
        pyth3 = False 
    else:
        raise ValueError

    apic = True

    if pyth3:
        import meshio
        mesh = meshio.read('/home/gary/temp/' + '3d__var_3_e20031120-070000-000.vtu')
        print(mesh)
        print(type(mesh.points))
        print(mesh.points.dtype)
        print(mesh.points.shape)

        p = mesh.point_data['p']
        #b1 = mesh.point_data['b1']
        #print(p)
        print(type(p))
        print(p.dtype)
        print(p.shape)
        #print(type(b1))
        #print(b1.dtype)
        #print(b1.shape)

        np.save("points_file.npy", mesh.points)

        if os.path.exists("B_file.npy"):
            mesh.point_data['b'][:,:] = np.load("B_file.npy", allow_pickle=apic)
            mesh.point_data['b1'][:,:] = np.load("B1_file.npy", allow_pickle=apic)
            mesh.point_data['j'][:,:] = np.load("JA_file.npy", allow_pickle=apic)
            print('writing mesh')
            mesh.write('3d__var_tmp.vtu')
            print('wrote mesh')

    else:
        points = np.load("points_file.npy")
        
        B_arr = B(points)
        B1_arr = B1(points)
        JA_arr = JA(points)

        np.save("B_file.npy", B_arr, allow_pickle=apic)
        np.save("B1_file.npy", B1_arr,allow_pickle=apic)
        np.save("JA_file.npy", JA_arr,allow_pickle=apic)
#create_files.py:96: RuntimeWarning: overflow encountered in cosh
#  ret[:, 1] = (-B_T/(delta*mu0))*np.cosh(X[:,2]/delta)**-2

