# todo test in python 3
import jitFORTRAN
import numpy as np


def demo1():
    n = 10
    c = np.zeros(n)
    a = np.ones(n)
    b = 2.*np.ones(n)

    # mind indentation, even though we are within function, we need to have 
    # the multiline string start at no tab to generate the correct string
    foo_fortran = '''!FOO
subroutine FOO(A,B,C,N)
    implicit none
    !will pass the transposed np array (which is in fortran order) having shape (3, N, M)

    integer, intent(in) :: N
    real(8), intent(in), dimension(0:N-1) :: A, B
!f2py intent(in) :: A, B
    real(8), intent(inout), dimension(0:N-1) :: C
!f2py intent(in,out) :: C

    integer :: i

    do i=0,N-1
        C(i) = C(i) + A(i) + B(i)
    end do
end subroutine FOO
'''

    foo_F = jitFORTRAN.Fortran_Subroutine(foo_fortran, 'FOO')

    print(c)
    foo_F.execute(a,b,c,n)
    print(c)
    foo_F.execute(a,b,c,n)
    print(c)

def demo2():
    #todo recreate large error of .../LUHMANN1979/create_files.py   by passing np.float32 to real(8)
    pass

def demo3():

    foo_fortran = '''!FOO
subroutine FOO(A,B,C,N)
    implicit none
    call 
end subroutine FOO
'''

    # check if f loop calling subroutine is slower than just floop (for mead)
    with open('mead.f') as f:
        lines = concat(f.readlines())#!!!

    foo_fortran = foo_fortran + lines

demo1()
demo2()
