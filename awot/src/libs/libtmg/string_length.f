      INTEGER FUNCTION STRING_LENGTH(STRING)

C  Thomas Matejka NOAA/NSSL 19 May 2000

C  This function returns the length of the string STRING to the last
C  non-blank character. If STRING contains only blanks, then the
C  function returns 0. If STRING contains no blanks, then the function
C  returns the length of STRING.

      IMPLICIT NONE
      CHARACTER(LEN=*)::STRING
      INTEGER::I,IEND

      IEND=LEN(STRING)
      DO I=IEND,1,-1
         IF(STRING(I:I).NE.' ')THEN
            STRING_LENGTH=I
            RETURN
         ENDIF
      ENDDO
      STRING_LENGTH=0

      END FUNCTION STRING_LENGTH
