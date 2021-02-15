import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../')
import jitFORTRAN

script = '''!MEADV
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
'''# does the !f2py actually do anything?

with open('mead.f', 'r') as include:
    script = script + ''.join(include.readlines())

print(script)

meadV_F = jitFORTRAN.Fortran_Subroutine(script, 'MEADV')

x = np.arange(10, dtype=np.float64)+1.
y = np.arange(10, dtype=np.float64)+1.
z = np.arange(10, dtype=np.float64)+1.

bx = np.nan*np.empty((10,), dtype=np.float64)
by = np.nan*np.empty((10,), dtype=np.float64)
bz = np.nan*np.empty((10,), dtype=np.float64)

meadV_F.execute(x, y, z, 1, bx, by, bz, 10)

print(bx, by, bz)
