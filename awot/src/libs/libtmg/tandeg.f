      REAL FUNCTION TANDEG(X)

C  Thomas Matejka NOAA/NSSL 23 Feburary 1993

C  This function returns the tangent of X, where X is in degrees.

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL::X

      TANDEG=TAN(X*RADDEG)

      END FUNCTION TANDEG
