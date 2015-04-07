      Function RADW(VLAT,VLON,U,V)

      R=SQRT(VLAT*VLAT + VLON*VLON)
      IF (R .eq. 0.0) Go To 1
      RADW=(VLON*U+VLAT*V)/R
      Return
    1 RADW=0.0
      Return
      End
