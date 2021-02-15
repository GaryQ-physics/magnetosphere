!SUBV
subroutine SUBV(x, N)
    implicit none
    integer(8), intent(in) :: N
!f2py intent(in) :: N
    real(4), intent(in), dimension(N) :: x
!f2py intent(in) :: X
    integer(8) :: ind

    do ind=1,N
        print*, x(ind)
! fortran 90
! print an array of numbers to screen
    end do

end subroutine SUBV
