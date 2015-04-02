      SUBROUTINE LOAD_FLOAT_TRUNC_RIGHT(A,NDECPT_MAX,OUTSTRING,NCHAR)

C  Thomas Matejka NOAA/NSSL 20 April 1993

C  This subroutine writes the floating-point number A right justified
C  into the string OUTSTRING.  The number is written with a maximum of
C  NDECPT_MAX digits after the decimal point.  Trailing zeroes are not
C  written.  The decimal point is not written if no digits follow it.
C  If OUTSTRING is longer than needed, the beginning of OUTSTRING is
C  filled with blanks.

C  The number of digits used to represent the number is returned as
C  NCHAR.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(MAX_STRING) FMT
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      LOGICAL NON_ZERO
      INTEGER NDECPT_MAX,NCHAR,I,N,IMIN,IMAX
      REAL A

      IF(NDECPT_MAX.LT.0)THEN
         N=0
      ELSE
         N=NDECPT_MAX
      ENDIF
      WRITE(FMT,*)'(F',MAX_NUMBER_STRING,'.',N,')'
      WRITE(STORE,FMT)A
      IF(STORE(1:1).NE.' ')THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'LOAD_FLOAT_TRUNC_RIGHT:  MEMORY ',
     $   'EXCEEDED.  INCREASE MAX_NUMBER_STRING.'
         STOP
      ENDIF

      NON_ZERO=.FALSE.
      DO 1 I=MAX_NUMBER_STRING,1,-1
         IF(STORE(I:I).EQ.'.')THEN
            IMAX=I-1
            GOTO 2
         ENDIF
         IF(STORE(I:I).NE.'0')THEN
            IMAX=I
            NON_ZERO=.TRUE.
            GOTO 2
         ENDIF
1     CONTINUE
2     CONTINUE

      DO 3 I=IMAX,1,-1
         IF(STORE(I:I).NE.' ')THEN
            GOTO 4
         ELSE
            NON_ZERO=.TRUE.
         ENDIF
3     CONTINUE
4     CONTINUE
      IMIN=I+1

      IF(NON_ZERO)THEN
         CALL LOAD_STRING_RIGHT(STORE(IMIN:IMAX),OUTSTRING)
         NCHAR=IMAX-IMIN+1
      ELSE
         CALL LOAD_STRING_RIGHT('0',OUTSTRING)
         NCHAR=1
      ENDIF
      RETURN
      END
