      SUBROUTINE LOAD_STRING_RIGHT(INSTRING,OUTSTRING)

C  Thomas Matejka NOAA/NSSL 26 March 1993

C  This subroutine copies the string INSTRING right justified into the
C  string OUTSTRING.  If OUTSTRING is longer than needed, the beginning
C  of OUTSTRING is filled with blanks.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER*(*) INSTRING,OUTSTRING
      INTEGER INEND,OUTEND,J,JMIN

      INEND=LEN(INSTRING)
      OUTEND=LEN(OUTSTRING)
      IF(OUTEND.LT.INEND)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'LOAD_STRING_RIGHT:  OUTSTRING ',
     $   'IS TOO SHORT.'
         STOP
      ENDIF
      JMIN=OUTEND-(INEND-1)
      OUTSTRING(JMIN:OUTEND)=INSTRING(1:INEND)
      IF(JMIN.GT.1)THEN
         DO 1 J=1,JMIN-1
            OUTSTRING(J:J)=' '
1        CONTINUE
      ENDIF
      RETURN
      END
