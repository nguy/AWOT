      REAL FUNCTION SINDEG(X)

C  Thomas Matejka NOAA/NSSL 23 February 1993

C  This function returns the sine of X, where X is in degrees.

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL::X

      SINDEG=SIN(X*RADDEG)

      END FUNCTION SINDEG
