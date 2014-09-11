      REAL FUNCTION AINTOUT(A)

C  Thomas Matejka NOAA/NSSL 18 February 1993

C  This function returns A rounded out to a whole number.

C  3.2 becomes 4., and -3.2 becomes -4..

      IMPLICIT NONE
      REAL A

      AINTOUT=AINT(A)
      IF(AINTOUT.NE.A)THEN
         AINTOUT=AINTOUT+SIGN(1.,A)
      ENDIF
      RETURN
      END
