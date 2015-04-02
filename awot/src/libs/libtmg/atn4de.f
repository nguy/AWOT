      REAL FUNCTION ATN4DE(A,B)

C  Thomas Matejka NOAA/NSSL 14 July 1995

C  This function returns the arctangent of A/B in degrees in the range
C  [0.,360.).

C  If A = 0. and B = 0., the arctangent is arbitrary and is returned 0.

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL A,B

      IF(A.NE.0..OR.
     $B.NE.0.)THEN
         ATN4DE=ATAN2(A,B)*DEGRAD
         IF(ATN4DE.LT.0.)THEN
            ATN4DE=ATN4DE+360.
         ENDIF
         IF(ATN4DE.GE.360.)THEN
            ATN4DE=0.
         ENDIF
      ELSE
         ATN4DE=0.
      ENDIF
      RETURN
      END
