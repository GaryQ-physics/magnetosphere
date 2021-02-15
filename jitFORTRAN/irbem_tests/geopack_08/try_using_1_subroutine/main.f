!##################### main.f ##########################################
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER*8, PARAMETER :: N = 5
      REAL*8 XGEO_V(N),YGEO_V(N),ZGEO_V(N),XMAG_V(N),YMAG_V(N),ZMAG_V(N)
      INTEGER*4 IYEAR,IDAY,IHOUR,MIN,ISEC, J

      IYEAR = 1
      IDAY = 1
      IHOUR = 1
      MIN = 1
      ISEC = 1

      J = 1

      XGEO_V = (/ 0., 1., 2., 3., 4. /)
      YGEO_V = (/ 0., 1., 2., 3., 4. /)
      ZGEO_V = (/ 0., 1., 2., 3., 4. /)
      XMAG_V = (/ 0., 1., 2., 3., 4. /)
      YMAG_V = (/ 0., 1., 2., 3., 4. /)
      ZMAG_V = (/ 0., 1., 2., 3., 4. /)

      print*,XGEO_V
      print*,YGEO_V
      print*,ZGEO_V
      print*,XMAG_V
      print*,YMAG_V
      print*,ZMAG_V

      CALL GEOMAG_08_V(XGEO_V,YGEO_V,ZGEO_V,XMAG_V,YMAG_V,ZMAG_V,
     1                       J,
     2                       IYEAR,IDAY,IHOUR,MIN,ISEC,
     3                       N)

      print*,XGEO_V
      print*,YGEO_V
      print*,ZGEO_V
      print*,XMAG_V
      print*,YMAG_V
      print*,ZMAG_V

      END
