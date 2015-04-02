      REAL FUNCTION C_FROM_K(K,BADDATA)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This function converts a temperature from Kelvin to Celcius.

C  The function returns BADDATA if K = BADDATA.

      IMPLICIT NONE
      REAL K,BADDATA

      IF(K.NE.BADDATA)THEN
         C_FROM_K=K-273.15
      ELSE
         C_FROM_K=BADDATA
      ENDIF
      RETURN
      END
