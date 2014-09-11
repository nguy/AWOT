      REAL FUNCTION K_FROM_F(F,BADDATA)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This function converts a temperature from Fahrenheit to Kelvin.

C  The function returns BADDATA if F = BADDATA.

      IMPLICIT NONE
      REAL F,BADDATA

      IF(F.NE.BADDATA)THEN
         K_FROM_F=(F-32.)*5./9.+273.15
      ELSE
         K_FROM_F=BADDATA
      ENDIF
      RETURN
      END
