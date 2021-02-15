
!#######################################################################
!#######   sub.f     ##########################################
!#######################################################################
! fortran 77
! print one imported float32 to screen 

      SUBROUTINE sub(x)
      IMPLICIT NONE
      REAL*4 x
      REAL*4 scaleby

      COMMON /blockname/scaleby

      print*, scaleby*x

      END

