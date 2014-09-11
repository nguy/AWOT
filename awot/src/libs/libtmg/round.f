      REAL FUNCTION ROUND(X,Y)

C  Thomas Matejka NOAA/NSSL 18 February 1993

C  This function returns X rounded to the nearest Y.

      IMPLICIT NONE
      REAL X,Y,YY

      IF(Y.EQ.0.)THEN
         ROUND=X
      ELSE
         YY=ABS(Y)
         ROUND=YY*AINT((X+YY*SIGN(0.5,X))/YY)
      ENDIF
      RETURN
      END
