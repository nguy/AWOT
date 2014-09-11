      SUBROUTINE LOAD_FLOAT(A,NDECPT,OUTSTRING,NCHAR)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine writes the floating-point number A left justified
C  into the string OUTSTRING.  The number is written with NDECPT digits
C  after the decimal point.  If OUTSTRING is longer than needed, the end
C  of OUTSTRING is filled with blanks.

C  The number of digits used to represent the number is returned as
C  NCHAR.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_STRING) FMT
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NDECPT,NCHAR,I,N,IMIN
      REAL A

      IF(NDECPT.LT.0)THEN
         N=0
      ELSE
         N=NDECPT
      ENDIF
      WRITE(FMT,*)'(F',MAX_NUMBER_STRING,'.',N,')'
      WRITE(STORE,FMT)A
      IF(STORE(1:1).NE.' ')THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'LOAD_FLOAT:  MEMORY EXCEEDED.  ',
     $   'INCREASE MAX_NUMBER_STRING.'
         STOP
      ENDIF

      DO 1 I=MAX_NUMBER_STRING,1,-1
         IF(STORE(I:I).EQ.' ')THEN
            GOTO 2
         ENDIF
1     CONTINUE
2     CONTINUE
      IMIN=I+1

      CALL LOAD_STRING(STORE(IMIN:MAX_NUMBER_STRING),OUTSTRING)
      NCHAR=MAX_NUMBER_STRING-IMIN+1
      RETURN
      END
