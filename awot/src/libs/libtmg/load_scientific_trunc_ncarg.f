      SUBROUTINE LOAD_SCIENTIFIC_TRUNC_NCARG(A,NDECPT_MAX,OUTSTRING,
     $NCHAR)

C  Thomas Matejka NOAA/NSSL 19 May 1994

C  This subroutine writes the floating-point number A left justified
C  into the string OUTSTRING.  The number is written in ncargraphics
C  scientific notation with a maximum of NDECPT_MAX digits after the
C  decimal point.  Trailing zeroes are not written.  The decimal point
C  is not written if no digits follow it.  If OUTSTRING is longer than
C  needed, the end of OUTSTRING is filled with blanks.

C  The number of digits used to represent the number is returned as
C  NCHAR.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER S_L
      CHARACTER*(MAX_NUMBER_STRING) STORE
      CHARACTER*(*) OUTSTRING
      INTEGER NDECPT_MAX,NCHAR,I,N,POS_POINT,POS_E
      REAL A,B

      IF(NDECPT_MAX.LT.0)THEN
         N=0
      ELSE
         N=NDECPT_MAX
      ENDIF
      B=0.1*A
      CALL LOAD_SCIENTIFIC_TRUNC(B,N+1,STORE,NCHAR)
      
      DO 1 I=MAX_NUMBER_STRING,1,-1
         IF(STORE(I:I).EQ.'.')THEN
            POS_POINT=I
            GOTO 2
         ENDIF
1     CONTINUE
2     CONTINUE
      STORE(POS_POINT:POS_POINT)=STORE(POS_POINT+1:POS_POINT+1)
      STORE(POS_POINT+1:POS_POINT+1)='.'
      DO 3 I=MAX_NUMBER_STRING,1,-1
         IF(STORE(I:I).EQ.'E')THEN
            POS_E=I
            GOTO 4
         ENDIF
3     CONTINUE
4     CONTINUE
      IF(POS_E.EQ.POS_POINT+2)THEN
         CALL LOAD_STRING(STORE(1:POS_E-2),OUTSTRING)
      ELSE
         CALL LOAD_STRING(STORE(1:POS_E-1),OUTSTRING)
      ENDIF
      CALL APPEND_STRING(0,':K: :P::0137::K: :P:10:S:',OUTSTRING)
      READ(STORE(POS_E+1:MAX_NUMBER_STRING),*)I
      CALL APPEND_INTEGER(0,I,OUTSTRING)
      CALL APPEND_STRING(0,':N:',OUTSTRING)
      NCHAR=S_L(OUTSTRING)

      RETURN
      END            
