
!#######################################################################
!#######   subroutine.f     ##########################################
!#######################################################################
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

! print one imported float32 to screen 
      SUBROUTINE sub(x)
      IMPLICIT NONE
      REAL*4 x
      REAL*4 scaleby

      COMMON /blockname/scaleby

      print*, scaleby*x

      END
