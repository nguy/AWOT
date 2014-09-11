      REAL FUNCTION ROUND_SIGNIF(A,NSIG)

C  Thomas Matejka NOAA/NSSL 18 February 1993

C  This function returns A rounded so as to include only NSIG
C  significant digits.

      IMPLICIT NONE
      REAL ROUND
      INTEGER NSIG,N,IEXP
      REAL A,B

      IF(NSIG.LE.0)THEN
         N=1
      ELSE
         N=NSIG
      ENDIF
      IF(A.NE.0.)THEN
         B=ALOG10(ABS(A))
         IEXP=IFIX(B)+1
         IF(B.LT.0..AND.
     $   AMOD(B,1.).NE.0.)THEN
            IEXP=IEXP-1
         ENDIF
      ELSE
         IEXP=1
      ENDIF
      ROUND_SIGNIF=ROUND(A/10.**IEXP,10.**(-N))*10.**IEXP
      RETURN
      END
