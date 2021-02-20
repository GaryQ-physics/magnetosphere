
!################ outer_sub.f #########################################

!MEADV
      SUBROUTINE MEADV(xxv, yyv, zzv, kp, abxv, abyv, abzv, N)
      IMPLICIT NONE
C inputed dimension of arrays
      INTEGER*8 N
C INPUT:
      REAL*8  tilt
      REAL*8  xxv(N), yyv(N), zzv(N)
!f2py intent(in) ::  xxv, yyv, zzv
      INTEGER*4 KP
C OUTPUT:
      REAL*8  abxv(N), abyv(N), abzv(N)
!f2py intent(in,out) :: abxv, abyv, abzv
C local:
      COMMON /dip_ang/tilt
      INTEGER*8 ind
C initialize common block tilt for use in MEAD subroutine:
      TILT = 1. ! even without setting, it doesnt crash for some reason

      DO 10 ind=1,N
        CALL MEAD(xxv(ind), yyv(ind), zzv(ind), kp, abxv(ind),
     1              abyv(ind), abzv(ind))
10    CONTINUE
      END
