      SUBROUTINE LOAD_SCIENTIFIC_RIGHT(A,NDECPT,OUTSTRING,NCHAR)

C  Thomas Matejka NOAA/NSSL 12 May 1994

C  This subroutine writes the floating-point number A right justified
C  into the string OUTSTRING.  The number is written in scientific
C  notation with NDECPT digits after the decimal point.  If OUTSTRING is
C  longer than needed, the beginning of OUTSTRING is filled with blanks.

C  The number of digits used to represent the number is returned as
C  NCHAR.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_STRING) FMT
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NDECPT,NCHAR,I,N,IMIN
      REAL A

      IF(NDECPT.LT.1)THEN
         N=1
      ELSE
         N=NDECPT
      ENDIF
      WRITE(FMT,*)'(E',MAX_NUMBER_STRING,'.',N,')'
      WRITE(STORE,FMT)A
      IF(STORE(1:1).NE.' ')THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'LOAD_SCIENTIFIC_RIGHT:  MEMORY ',
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
