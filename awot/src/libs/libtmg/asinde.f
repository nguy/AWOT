      REAL FUNCTION ASINDE(XIN)

C  Thomas Matejka NOAA/NSSL 14 July 1995

C  This function returns the arcsine of XIN in degrees in the range
C  [-90.,90.].

C  If XIN is greater than 1. or less than -1., then XIN is assumed to be
C  1.  or -1..

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL XIN,X

      IF(ABS(XIN).GT.1.)THEN
         X=SIGN(1.,XIN)
      ELSE
         X=XIN
      ENDIF
      ASINDE=ASIN(X)*DEGRAD
      RETURN
      END
