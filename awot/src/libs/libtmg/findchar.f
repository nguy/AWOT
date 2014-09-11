      SUBROUTINE FINDCHAR(STRING,CHAR,OCCURRENCE,POS,SUCCESS)

C  Thomas Matejka NOAA/NSSL 19 May 2000

C  This subroutine calculates the position of the OCCURRENCEth
C  occurrence of the character CHAR in the string STRING.

C  SUCCESS is returned .FALSE. if and only if OCCURRENCE occurrences of
C  CHAR cannot be found.

      IMPLICIT NONE
      CHARACTER::CHAR
      CHARACTER(LEN=*)::STRING
      LOGICAL::SUCCESS
      INTEGER::I,IEND,OCCURRENCE,FOUND,POS

      IEND=LEN(STRING)
      FOUND=0
      DO I=1,IEND
         IF(STRING(I:I).EQ.CHAR)THEN
            FOUND=FOUND+1
            IF(FOUND.EQ.OCCURRENCE)THEN
               SUCCESS=.TRUE.
               POS=I
               RETURN
            ENDIF
         ENDIF
      ENDDO
      SUCCESS=.FALSE.

      END SUBROUTINE FINDCHAR
