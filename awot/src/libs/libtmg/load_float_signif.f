      SUBROUTINE LOAD_FLOAT_SIGNIF(A,NSIG,OUTSTRING,NCHAR)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine writes the floating-point number A left justified
C  into the string OUTSTRING.  The number is written with NSIG
C  significant digits.  If OUTSTRING is longer than needed, the end of
C  OUTSTRING is filled with blanks.

C  The number of digits used to represent the number is returned as
C  NCHAR.

      IMPLICIT NONE
      REAL ROUND
      CHARACTER*(*) OUTSTRING
      INTEGER NSIG,NCHAR,IEXP,N,NDECPT
      REAL A,B,C

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
      C=A/10.**IEXP
      C=ROUND(C,10.**(-N))
      C=C*10.**IEXP
      IF(N-IEXP.GE.0)THEN
         NDECPT=N-IEXP
      ELSE
         NDECPT=0
      ENDIF
      CALL LOAD_FLOAT(C,NDECPT,OUTSTRING,NCHAR)
      RETURN
      END
