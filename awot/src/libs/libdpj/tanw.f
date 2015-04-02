      Function TANW(VLAT,VLON,U,V)

      R=SQRT(VLAT*VLAT+VLON*VLON)
      IF (R .eq. 0.0) Go To 1
      TANW=(VLON*V-VLAT*U)/R
      Return
    1 TANW=0.0
      Return
      End
