      REAL FUNCTION K_FROM_C(C,BADDATA)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This function converts a temperature from Celcius to Kelvin.

C  The function returns BADDATA if C = BADDATA.

      IMPLICIT NONE
      REAL C,BADDATA

      IF(C.NE.BADDATA)THEN
         K_FROM_C=C+273.15
      ELSE
         K_FROM_C=BADDATA
      ENDIF
      RETURN
      END
