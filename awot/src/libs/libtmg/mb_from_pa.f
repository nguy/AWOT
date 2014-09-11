      REAL FUNCTION MB_FROM_PA(PA,BADDATA)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This function converts a pressure from Pascals to millibars.

C  The function returns BADDATA if PA = BADDATA.

      IMPLICIT NONE
      REAL PA,BADDATA

      IF(PA.NE.BADDATA)THEN
         MB_FROM_PA=PA*0.01
      ELSE
         MB_FROM_PA=BADDATA
      ENDIF
      RETURN
      END
