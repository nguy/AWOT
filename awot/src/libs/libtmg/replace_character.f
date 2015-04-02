      SUBROUTINE REPLACE_CHARACTER(STRING,OLD_CHAR,NEW_CHAR)

C  Thomas Matejka NOAA/NSSL 15 March 1993

C  This subroutine replaces all occurrences of the character OLD_CHAR
C  with the character NEW_CHAR in the string STRING.

      IMPLICIT NONE
      CHARACTER OLD_CHAR,NEW_CHAR
      CHARACTER*(*) STRING
      INTEGER IEND,I

      IEND=LEN(STRING)
      DO 1 I=1,IEND
         IF(STRING(I:I).EQ.OLD_CHAR)THEN
            STRING(I:I)=NEW_CHAR
         ENDIF
1     CONTINUE
      RETURN
      END
