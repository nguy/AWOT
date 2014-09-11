      REAL FUNCTION F_FROM_K(K,BADDATA)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This function converts a temperature from Kelvin to Fahrenheit.

C  The function returns BADDATA if K = BADDATA.

      IMPLICIT NONE
      REAL K,BADDATA

      IF(K.NE.BADDATA)THEN
         F_FROM_K=(K-273.15)*9./5.+32.
      ELSE
         F_FROM_K=BADDATA
      ENDIF
      RETURN
      END
