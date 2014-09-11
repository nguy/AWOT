      INTEGER FUNCTION S_L(STRING)

C  Thomas Matejka NOAA/NSSL 7 June 2002

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
            S_L=I
            RETURN
         ENDIF
      ENDDO
      S_L=0

      END FUNCTION S_L
