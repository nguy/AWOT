      SUBROUTINE APPEND_FLOAT_TRUNC(NBLANKS,A,NDECPT_MAX,OUTSTRING)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine writes NBLANKS blanks and the floating-point number A
C  left justified into the string OUTSTRING starting after the last
C  non-blank character.  The number is written with a maximum of
C  NDECPT_MAX digits after the decimal point.  Trailing zeroes are not
C  written.  The decimal point is not written if no digits follow it.
C  The rest of OUTSTRING is filled with blanks.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NBLANKS,NDECPT_MAX,NCHAR
      REAL A

      CALL LOAD_FLOAT_TRUNC(A,NDECPT_MAX,STORE,NCHAR)
      CALL APPEND_STRING(NBLANKS,STORE(1:NCHAR),OUTSTRING)
      RETURN
      END
