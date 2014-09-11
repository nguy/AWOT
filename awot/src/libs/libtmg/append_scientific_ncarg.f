      SUBROUTINE APPEND_SCIENTIFIC_NCARG(NBLANKS,A,NDECPT,OUTSTRING)

C  Thomas Matejka NOAA/NSSL 13 May 1994

C  This subroutine writes NBLANKS blanks and the floating-point number A
C  left justified into the string OUTSTRING starting after the last
C  non-blank character.  The number is written in ncargraphics
C  scientific notation with NDECPT digits after the decimal point.  The
C  rest of OUTSTRING is filled with blanks.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NBLANKS,NDECPT,NCHAR
      REAL A

      CALL LOAD_SCIENTIFIC_NCARG(A,NDECPT,STORE,NCHAR)
      CALL APPEND_STRING(NBLANKS,STORE(1:NCHAR),OUTSTRING)
      RETURN
      END
