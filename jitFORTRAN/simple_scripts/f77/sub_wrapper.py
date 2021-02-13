import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../')
import jitFORTRAN

script = '''!SUBV
      SUBROUTINE subV(x,N)
      IMPLICIT NONE
      INTEGER*8 N
      REAL*4 x(N)
!f2py intent(in) :: x
      INTEGER*8 ind

      DO 10 ind=1,N
        call sub(x(ind))
10    CONTINUE
      END
'''# does the !f2py actually do anything?

with open('sub.f', 'r') as include:
    script = script + ''.join(include.readlines())

print(script)

subV_F = jitFORTRAN.Fortran_Subroutine(script, 'SUBV')

x = np.arange(10, dtype=np.float32)
subV_F.execute(x)
a = np.arange(10, dtype=np.float32)**2
subV_F.execute(a, 10)
