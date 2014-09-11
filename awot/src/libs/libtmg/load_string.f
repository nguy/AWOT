      SUBROUTINE LOAD_STRING(INSTRING,OUTSTRING)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine copies the string INSTRING left justified into the
C  string OUTSTRING.  If OUTSTRING is longer than needed, the end of
C  OUTSTRING is filled with blanks.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(*) INSTRING,OUTSTRING
      INTEGER INEND,OUTEND,J,JMAX

      INEND=LEN(INSTRING)
      OUTEND=LEN(OUTSTRING)
      IF(OUTEND.LT.INEND)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'LOAD_STRING:  OUTSTRING IS TOO ',
     $   'SHORT.'
         STOP
      ENDIF
      JMAX=INEND
      OUTSTRING(1:JMAX)=INSTRING(1:INEND)
      IF(JMAX.LT.OUTEND)THEN
         DO 1 J=JMAX+1,OUTEND
            OUTSTRING(J:J)=' '
1        CONTINUE
      ENDIF
      RETURN
      END
