
!################ outer_sub.f #########################################

!GEOMAG_08_V
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
