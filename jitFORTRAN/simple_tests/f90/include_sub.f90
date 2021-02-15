! fortran 90
! print an array of numbers to screen by 
!   repeadedly calling sub by
!     including sub.f90 source code
!https://stackoverflow.com/questions/28134327/fortran-split-module-into-multiple-files
!https://stackoverflow.com/questions/45223705/how-to-compile-fortran-function-subroutine-in-separate-files-into-a-single-modul

subroutine subV(x, N)
    implicit none
    integer(8) :: N
    real(4), dimension(N) :: x
!f2py intent(in) :: x
    integer(8) :: ind

    do ind=1,N
        call sub(x(ind))
    end do

  contains 
    include 'sub.f90'

end subroutine subV

