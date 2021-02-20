
!##################### main.f ##########################################
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER*8, PARAMETER :: N = 5
      REAL*8 X_V(N),Y_V(N),Z_V(N),abxv(N),abyv(N),abzv(N)
      INTEGER*4 KP

      KP = 1

      X_V =  (/ 0., 1., 2., 3., 4. /)
      Y_V =  (/ 0., 1., 2., 3., 4. /)
      Z_V =  (/ 0., 1., 2., 3., 4. /)
      abxv = (/ 0., 1., 2., 3., 4. /)
      abxv = (/ 0., 1., 2., 3., 4. /)
      abxv = (/ 0., 1., 2., 3., 4. /)

      print*,X_V
      print*,Y_V
      print*,Z_V
      print*,abxv
      print*,abyv
      print*,abzv

      CALL MEADV(X_V, Y_V, Z_V, kp, abxv, abyv, abzv, N)

      print*,X_V
      print*,Y_V
      print*,Z_V
      print*,abxv
      print*,abyv
      print*,abzv

      END
