      REAL FUNCTION DBX_TO_X(DBX)

C  Thomas Matejka NOAA/NSSL 18 February 1993

C  This function converts DBX from decibels to a linear value.

      IMPLICIT NONE
      REAL DBX

      DBX_TO_X=10.**(DBX/10.)
      RETURN
      END
