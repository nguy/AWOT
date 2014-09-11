      SUBROUTINE APPEND_SCIENTIFIC_TRUNC(NBLANKS,A,NDECPT_MAX,OUTSTRING)

C  Thomas Matejka NOAA/NSSL 12 May 1994

C  This subroutine writes NBLANKS blanks and the floating-point number A
C  left justified into the string OUTSTRING starting after the last
C  non-blank character.  The number is written in scientific notation
C  with a maximum of NDECPT_MAX digits after the decimal point.
C  Trailing zeroes are not written.  The rest of OUTSTRING is filled
C  with blanks.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NBLANKS,NDECPT_MAX,NCHAR
      REAL A

      CALL LOAD_SCIENTIFIC_TRUNC(A,NDECPT_MAX,STORE,NCHAR)
      CALL APPEND_STRING(NBLANKS,STORE(1:NCHAR),OUTSTRING)
      RETURN
      END
