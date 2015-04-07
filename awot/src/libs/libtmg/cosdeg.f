      REAL FUNCTION COSDEG(X)

C  Thomas Matejka NOAA/NSSL 23 February 1993

C  This function returns the cosine of X, where X is in degrees.

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL::X

      COSDEG=COS(X*RADDEG)

      END FUNCTION COSDEG
