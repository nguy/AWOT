      SUBROUTINE EXCISE_STRING(STRING,ISTART,IEND)

C  Thomas Matejka NOAA/NSSL 16 October 1996

C  This subroutine excises character numbers ISTART to IEND in the
C  character string STRING.  Trailing blanks are added to fill STRING.

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER(LEN=*)::STRING
      INTEGER::ISTART,IEND,L

      L=LEN(STRING)
      IF(IEND.LT.L)THEN
         STRING(ISTART:L-IEND+ISTART-1)=STRING(IEND+1:L)
      ENDIF
      STRING(L-IEND+ISTART:L)=''
      END


