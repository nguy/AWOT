      SUBROUTINE APPEND_SCIENTIFIC(NBLANKS,A,NDECPT,OUTSTRING)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine writes NBLANKS blanks and the floating-point number A
C  left justified into the string OUTSTRING starting after the last
C  non-blank character.  The number is written in scientific notation
C  with NDECPT digits after the decimal point.  The rest of OUTSTRING is
C  filled with blanks.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NBLANKS,NDECPT,NCHAR
      REAL A

      CALL LOAD_SCIENTIFIC(A,NDECPT,STORE,NCHAR)
      CALL APPEND_STRING(NBLANKS,STORE(1:NCHAR),OUTSTRING)
      RETURN
      END
