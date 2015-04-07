      SUBROUTINE LOAD_SCIENTIFIC_TRUNC(A,NDECPT_MAX,OUTSTRING,NCHAR)

C  Thomas Matejka NOAA/NSSL 12 May 1994

C  This subroutine writes the floating-point number A left justified
C  into the string OUTSTRING.  The number is written in scientific
C  notation with a maximum of NDECPT_MAX digits after the decimal point.
C  Trailing zeroes are not written.  If OUTSTRING is longer than needed,
C  the end of OUTSTRING is filled with blanks.

C  The number of digits used to represent the number is returned as
C  NCHAR.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_STRING) FMT
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NDECPT_MAX,NCHAR,I,N,IMIN,POS_E,POS_POINT
      REAL A

      IF(NDECPT_MAX.LT.1)THEN
         N=1
      ELSE
         N=NDECPT_MAX
      ENDIF
      WRITE(FMT,*)'(E',MAX_NUMBER_STRING,'.',N,')'
      WRITE(STORE,FMT)A
      IF(STORE(1:1).NE.' ')THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'LOAD_SCIENTIFIC:  MEMORY ',
     $   'EXCEEDED.  INCREASE MAX_NUMBER_STRING.'
         STOP
      ENDIF

      DO 1 I=MAX_NUMBER_STRING,1,-1
         IF(STORE(I:I).EQ.'E')THEN
            POS_E=I
            GOTO 2
         ENDIF
1     CONTINUE
2     CONTINUE
      DO 3 I=POS_E-1,1,-1
         IF(STORE(I:I).EQ.'.')THEN
            POS_POINT=I
            GOTO 4
         ENDIF
3     CONTINUE
4     CONTINUE
      DOWHILE(POS_E-1.GE.POS_POINT+2)
         IF(STORE(POS_E-1:POS_E-1).EQ.'0')THEN
            DO 5 I=POS_E-1,2,-1
               STORE(I:I)=STORE(I-1:I-1)
5           CONTINUE
            STORE(1:1)=' '
            POS_POINT=POS_POINT+1
         ELSE
            GOTO 6
         ENDIF
      ENDDO
6     CONTINUE

      DO 7 I=MAX_NUMBER_STRING,1,-1
         IF(STORE(I:I).EQ.' ')THEN
            GOTO 8
         ENDIF
7     CONTINUE
8     CONTINUE
      IMIN=I+1

      CALL LOAD_STRING(STORE(IMIN:MAX_NUMBER_STRING),OUTSTRING)
      NCHAR=MAX_NUMBER_STRING-IMIN+1
      RETURN
      END







