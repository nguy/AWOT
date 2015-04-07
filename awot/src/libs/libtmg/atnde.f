      REAL FUNCTION ATNDE(A,B)

C  Thomas Matejka NOAA/NSSL 14 July 1995

C  This function returns the arctangent of A/B in degrees in the range
C  [-90.,90.].  The function assumes that B is supposed to be
C  non-negative, even though it may be negative because of computational
C  error.  Therefore, arctangents less than -90. are set to -90., and
C  arctangents greater than 90. are set to 90..

C  If A = 0. and B = 0., the arctangent is arbitrary and is returned 0.

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL A,B

      IF(A.NE.0..OR.
     $B.NE.0.)THEN
         ATNDE=ATAN2(A,B)*DEGRAD
         IF(ATNDE.LT.-90.)THEN
            ATNDE=-90.
         ENDIF
         IF(ATNDE.GT.90.)THEN
            ATNDE=90.
         ENDIF
      ELSE
         ATNDE=0.
      ENDIF
      RETURN
      END
