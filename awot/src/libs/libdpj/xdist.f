      Function XDIST(XLON,CLON,PLAT)

C*** Returns DISTANCE (+/-) OF THE POINT(XLON,PLAT)  FROM THE ORIGIN
C*** (CLON,PLAT) ON THE CONSTANT LATITUDE LINE ,PLAT

      CONS = 0.017453292
      RCOS = 6378.388 * Cos(PLAT*CONS)
      XM = 1.0
      If (XLON.lt.CLON) XM = -1.0
      DL  =    ABS( CLON - XLON ) * CONS
      XDIST   = RCOS * DL *XM
      Return
      End
