      REAL FUNCTION INHG_FROM_PA(PA,BADDATA)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This function converts a pressure from Pascals to inches of mercury.

C  The function returns BADDATA if PA = BADDATA.

      IMPLICIT NONE
      REAL PA,BADDATA

      IF(PA.NE.BADDATA)THEN
         INHG_FROM_PA=PA/3386.39
      ELSE
         INHG_FROM_PA=BADDATA
      ENDIF
      RETURN
      END
