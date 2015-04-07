      REAL FUNCTION ATN3DE(A,B)

C  Thomas Matejka NOAA/NSSL 14 July 1995

C  This function returns the arctangent of A/B in degrees in the range
C  [-180.,180.).

C  If A = 0. and B = 0., the arctangent is arbitrary and is returned 0.

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL A,B

      IF(A.NE.0..OR.
     $B.NE.0.)THEN
         ATN3DE=ATAN2(A,B)*DEGRAD
         IF(ATN3DE.LT.-180.)THEN
            ATN3DE=ATN3DE+360.
         ENDIF
         IF(ATN3DE.GE.180.)THEN
            ATN3DE=-180.
         ENDIF
      ELSE
         ATN3DE=0.
      ENDIF
      RETURN
      END
