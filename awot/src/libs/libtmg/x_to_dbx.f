      REAL FUNCTION X_TO_DBX(X)

C  Thomas Matejka NOAA/NSSL 18 Februrary 1993

C  This function converts X from a linear value to decibels.

      IMPLICIT NONE
      REAL X

      X_TO_DBX=10.*ALOG10(X)
      RETURN
      END
