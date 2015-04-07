      SUBROUTINE LOAD_INTEGER_RIGHT(K,OUTSTRING,NCHAR)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine writes the integer K right justified into the string
C  OUTSTRING.  If OUTSTRING is longer than needed, the beginning of
C  OUTSTRING is filled with blanks.

C  The number of digits used to represent the number is returned as
C  NCHAR.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_STRING) FMT
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER K,NCHAR,I,IMIN

      WRITE(FMT,*)'(I',MAX_NUMBER_STRING,')'
      WRITE(STORE,FMT)K
      IF(STORE(1:1).NE.' ')THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'LOAD_INTEGER_RIGHT:  MEMORY ',
     $   'EXCEEDED.  INCREASE MAX_NUMBER_STRING.'
         STOP
      ENDIF

      DO 1 I=MAX_NUMBER_STRING,1,-1
         IF(STORE(I:I).EQ.' ')THEN
            GOTO 2
         ENDIF
1     CONTINUE
2     CONTINUE
      IMIN=I+1

      CALL LOAD_STRING_RIGHT(STORE(IMIN:MAX_NUMBER_STRING),OUTSTRING)
      NCHAR=MAX_NUMBER_STRING-IMIN+1
      RETURN
      END
