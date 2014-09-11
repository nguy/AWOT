      SUBROUTINE APPEND_FLOAT_SIGNIF_TRUNC(NBLANKS,A,NSIG_MAX,OUTSTRING)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine writes NBLANKS blanks and the floating-point number A
C  left justified into the string OUTSTRING starting after the last
C  non-blank character.  The number is written with a maximum of
C  NSIG_MAX significant digits.  Trailing zeroes are not written.  The
C  decimal point is not written if no digits follow it.  The rest of
C  OUTSTRING is filled with blanks.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NBLANKS,NSIG_MAX,NCHAR
      REAL A

      CALL LOAD_FLOAT_SIGNIF_TRUNC(A,NSIG_MAX,STORE,NCHAR)
      CALL APPEND_STRING(NBLANKS,STORE(1:NCHAR),OUTSTRING)
      RETURN
      END
