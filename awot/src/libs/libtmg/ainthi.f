      REAL FUNCTION AINTHI(A)

C  Thomas Matejka NOAA/NSSL 18 February 1993

C  This function returns A rounded up to a whole number.

C  3.2 becomes 4., and -3.2 becomes -3..

      IMPLICIT NONE
      REAL A

      AINTHI=AINT(A)
      IF(A.GT.0..AND.
     $AINTHI.NE.A)THEN
         AINTHI=AINTHI+1.
      ENDIF
      RETURN
      END
