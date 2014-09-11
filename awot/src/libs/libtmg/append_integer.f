      SUBROUTINE APPEND_INTEGER(NBLANKS,K,OUTSTRING)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine writes NBLANKS blanks and the integer K left
C  justified into the string OUTSTRING starting after the last non-blank
C  character.  The rest of OUTSTRING is filled with blanks.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NBLANKS,K,NCHAR

      CALL LOAD_INTEGER(K,STORE,NCHAR)
      CALL APPEND_STRING(NBLANKS,STORE(1:NCHAR),OUTSTRING)
      RETURN
      END
