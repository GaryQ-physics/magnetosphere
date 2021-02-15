!! fortran 90
!! main program for running the .f90 subroutines

!program main
!    implicit none
!    real(4) :: x = 5.
!
!    call sub(x)
!
!  contains
!    include 'sub.f90' 
!
!end program

!!!!!!!!!!!! this one does work, using just $ gfortran main.f90
!program main
!    implicit none
!    integer(8), parameter :: N=5
!    real(4), dimension(N) :: x = (/ 0., 1., 2., 3., 4. /)
!
!    call subV(x, N)
!
!  contains
!    include 'sub_arr.f90' 
!
!end program

!!!!!!!!!!!! this one  DOESNT WORK, due to multiple subroutines keywords clashing
!program main
!    implicit none
!    integer(8), parameter :: N=5
!    real(4), dimension(N) :: x = (/ 0., 1., 2., 3., 4. /)
!
!    call subV(x, N)
!
!  contains
!    include 'include_sub.f90' 
!
!end program


!!!!!!!!!!!!! build as either (both work):
!!!!!!!!!!!!!         $ gfortran main.f90 sub_arr.f90
!!!!!!!!!!!!!         $ gfortran main.f90 include_sub.f90
program main
    implicit none
    integer(8), parameter :: N=5
    real(4), dimension(N) :: x = (/ 0., 1., 2., 3., 4. /)

    call subV(x, N)

end program


