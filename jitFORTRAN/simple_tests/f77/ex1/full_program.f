
!#######################################################################
!#######   full_program.f     ##########################################
!#######################################################################

      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER*8, PARAMETER :: N = 5
      REAL*4 X_V(N),Y_V(N),Z_V(N)

      X_V = (/ 0., 1., 2., 3., 4. /)

      CALL subV(X_V,
     1            N
     2              )

      END
!SUBV
      SUBROUTINE subV(x,N)
      IMPLICIT NONE
      INTEGER*8 N
      REAL*4 x(N)
!f2py intent(in) :: x
      INTEGER*8 ind
      REAL*4 scaleby

      COMMON /blockname/scaleby
      scaleby = 2.

      DO 10 ind=1,N
        call sub(x(ind))
10    CONTINUE
      END
! fortran 77
! print one imported float32 to screen 

      SUBROUTINE sub(x)
      IMPLICIT NONE
      REAL*4 x
      REAL*4 scaleby

      COMMON /blockname/scaleby

      print*, scaleby*x

      END

