      SUBROUTINE LOAD_FLOAT_SIGNIF_TRUNC_RIGHT(A,NSIG_MAX,OUTSTRING,
     $NCHAR)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine writes the floating-point number A right justified
C  into the string OUTSTRING.  The number is written with a maximum of
C  NSIG_MAX significant digits.  Trailing zeroes are not written.  The
C  decimal point is not written if no digits follow it.  If OUTSTRING is
C  longer than needed, the beginning of OUTSTRING is filled with blanks.

C  The number of digits used to represent the number is returned as
C  NCHAR.

      IMPLICIT NONE
      REAL ROUND
      CHARACTER*(*) OUTSTRING
      INTEGER NSIG_MAX,NCHAR,IEXP,N,NDECPT_MAX
      REAL A,B,C

      IF(NSIG_MAX.LE.0)THEN
         N=1
      ELSE
         N=NSIG_MAX
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
         NDECPT_MAX=N-IEXP
      ELSE
         NDECPT_MAX=0
      ENDIF
      CALL LOAD_FLOAT_TRUNC_RIGHT(C,NDECPT_MAX,OUTSTRING,NCHAR)
      RETURN
      END
