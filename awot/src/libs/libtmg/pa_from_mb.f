      REAL FUNCTION PA_FROM_MB(MB,BADDATA)

C  Thomas Matejka NOAA/NSSL 11 March 1994

C  This function converts a pressure from millibars to Pascals.

C  The function returns BADDATA if MB = BADDATA.

      IMPLICIT NONE
      REAL MB,BADDATA

      IF(MB.NE.BADDATA)THEN
         PA_FROM_MB=MB*100.
      ELSE
         PA_FROM_MB=BADDATA
      ENDIF
      RETURN
      END
