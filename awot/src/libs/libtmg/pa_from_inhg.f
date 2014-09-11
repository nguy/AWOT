      REAL FUNCTION PA_FROM_INHG(INHG,BADDATA)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This function converts a pressure from inches of mercury to Pascals.

C  The function returns BADDATA if INHG = BADDATA.

      IMPLICIT NONE
      REAL INHG,BADDATA

      IF(INHG.NE.BADDATA)THEN
         PA_FROM_INHG=INHG*3386.39
      ELSE
         PA_FROM_INHG=BADDATA
      ENDIF
      RETURN
      END
