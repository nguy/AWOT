      SUBROUTINE APPEND_SCIENTIFIC_TRUNC_NCARG(NBLANKS,A,NDECPT_MAX,
     $OUTSTRING)

C  Thomas Matejka NOAA/NSSL 13 May 1994

C  This subroutine writes NBLANKS blanks and the floating-point number A
C  left justified into the string OUTSTRING starting after the last
C  non-blank character.  The number is written in ncargraphics
C  scientific notation with a maximum of NDECPT_MAX digits after the
C  decimal point.  Trailing zeroes are not written.  The rest of
C  OUTSTRING is filled with blanks.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NBLANKS,NDECPT_MAX,NCHAR
      REAL A

      CALL LOAD_SCIENTIFIC_TRUNC_NCARG(A,NDECPT_MAX,STORE,NCHAR)
      CALL APPEND_STRING(NBLANKS,STORE(1:NCHAR),OUTSTRING)
      RETURN
      END
