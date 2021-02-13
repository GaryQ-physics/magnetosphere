!#######################################################################
!############ this file may not compile all at once ####################
!#######################################################################


! subroutine.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SUBROUTINE subV(x,N)
      IMPLICIT NONE
      INTEGER*8 N
      REAL*4 x(N)
      INTEGER*8 ind

      DO 10 ind=1,N
        print*, x(ind)
10    CONTINUE

      END
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! main.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PROGRAM main
      IMPLICIT NONE
C     is this really fortran 77 ?
      INTEGER*8, PARAMETER :: N = 5
C     PARAMETER (N=5) # https://web.stanford.edu/class/me200c/tutorial_77/05_variables.html
      REAL*4 x(N)
      x = (/ 0., 1., 2., 3., 4. /)
  
      CALL subV(x, N)
      END
 
      INCLUDE 'subroutine.f'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!$ gfortran main.f
!$ ./a.out


!#######################################################################
!#######################################################################


! subroutine.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SUBROUTINE subV(x,N)
      IMPLICIT NONE
      INTEGER*8 N
      REAL*4 x(N)
      INTEGER*8 ind

      DO 10 ind=1,N
        print*, x(ind)
10    CONTINUE

      END
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! main.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PROGRAM main
      IMPLICIT NONE
      INTEGER*8, PARAMETER :: N = 5
      REAL*4 x(N)
      x = (/ 0., 1., 2., 3., 4. /)
  
      CALL subV(x, N)
      END
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!$ gfortran main.f subroutine.f
!$ ./a.out


!#######################################################################
!#######################################################################


! subroutine.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SUBROUTINE subV(x,N)
      IMPLICIT NONE
      INTEGER*8 N
      REAL*4 x(N)
      INTEGER*8 ind

      DO 10 ind=1,N
        call sub(x(ind))
10    CONTINUE
      END
      INCLUDE 'sub.f'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! main.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PROGRAM main
      IMPLICIT NONE
      INTEGER*8, PARAMETER :: N = 5
      REAL*4 x(N)
      x = (/ 0., 1., 2., 3., 4. /)
  
      CALL subV(x, N)
      END
 
      INCLUDE 'subroutine.f'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!$ gfortran main.f
!$ ./a.out


!#######################################################################
!#######################################################################


! subroutine.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SUBROUTINE subV(x,N)
      IMPLICIT NONE
      INTEGER*8 N
      REAL*4 x(N)
      INTEGER*8 ind

      DO 10 ind=1,N
        call sub(x(ind))
10    CONTINUE
      END
      INCLUDE 'sub.f'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! main.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      PROGRAM main
      IMPLICIT NONE
      INTEGER*8, PARAMETER :: N = 5
      REAL*4 x(N)
      x = (/ 0., 1., 2., 3., 4. /)
  
      CALL subV(x, N)
      END
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!$ gfortran main.f subroutine.f
!$ ./a.out
