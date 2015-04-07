      REAL FUNCTION ACOSDE(XIN)

C  Thomas Matejka NOAA/NSSL 14 July 1995

C  This function returns the arccosine of XIN in degrees in the range
C  [0.,180.].

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
      ACOSDE=ACOS(X)*DEGRAD
      RETURN
      END
