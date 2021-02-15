import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../')
import jitFORTRAN


script = '''!GEOMAG_08_V
      SUBROUTINE GEOMAG_08_V(XGEO_V,YGEO_V,ZGEO_V,XMAG_V,YMAG_V,ZMAG_V,
     1                       J,
     2                       IYEAR,IDAY,IHOUR,MIN,ISEC,
     3                       N)
C          vectorizes GEOMAG_08

C-----time of coordinate transform:
C     IYEAR   -  YEAR NUMBER (FOUR DIGITS)
C     IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
C     IHOUR -  HOUR OF DAY (00 TO 23)
C     MIN   -  MINUTE OF HOUR (00 TO 59)
C     ISEC  -  SECONDS OF MINUTE (00 TO 59)

C                    J>0                       J<0
C-----INPUT:  J,XGEO,YGEO,ZGEO           J,XMAG,YMAG,ZMAG
C-----OUTPUT:    XMAG,YMAG,ZMAG           XGEO,YGEO,ZGEO

      IMPLICIT NONE
C inputed dimension of arrays
      INTEGER*8 N
      REAL*8  XGEO_V(N), YGEO_V(N), ZGEO_V(N)
!f2py intent(in,out) ::  XGEO_V, YGEO_V, ZGEO_V
      REAL*8  XMAG_V(N), YMAG_V(N), ZMAG_V(N)
!f2py intent(in,out) ::  XMAG_V, YMAG_V, ZMAG_V
      INTEGER*4 J, IYEAR,IDAY,IHOUR,MIN,ISEC
      INTEGER*8 ind

C  RECALC_08 prepares elements of rot matrix and puts in common block
C  note since GEO and MAG are independent of solar wind, set 
C   can set VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0
      CALL RECALC_08(IYEAR,IDAY,IHOUR,MIN,ISEC,-4.0d2,0.0d0,0.0d0)

      DO 10 ind=1,N
      CALL GEOMAG_08(XGEO_V(ind),YGEO_V(ind),ZGEO_V(ind),
     1               XMAG_V(ind),YMAG_V(ind),ZMAG_V(ind), J)
10    CONTINUE
      END
'''

with open('geopack_08.f', 'r') as include:
    script = script + ''.join(include.readlines())

#print(script)

GEOMAG_08_V_F = jitFORTRAN.Fortran_Subroutine(script, 'GEOMAG_08_V')
GEOMAG_08_V_F.compile()
#x = np.arange(10, dtype=np.float64)+1.
#y = np.arange(10, dtype=np.float64)+1.
#z = np.arange(10, dtype=np.float64)+1.
#
#bx = np.nan*np.empty((10,), dtype=np.float64)
#by = np.nan*np.empty((10,), dtype=np.float64)
#bz = np.nan*np.empty((10,), dtype=np.float64)
#
#meadV_F.execute(x, y, z, 1, bx, by, bz, 10)
#
#print(bx, by, bz)
