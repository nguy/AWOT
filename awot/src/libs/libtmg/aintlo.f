      REAL FUNCTION AINTLO(A)

C  Thomas Matejka NOAA/NSSL 18 February 1993

C  This function returns A rounded down to a whole number.

C  3.2 becomes 3., and -3.2 becomes -4..

      IMPLICIT NONE
      REAL A

      AINTLO=AINT(A)
      IF(A.LT.0..AND.
     $AINTLO.NE.A)THEN
         AINTLO=AINTLO-1.
      ENDIF
      RETURN
      END
