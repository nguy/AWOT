      SUBROUTINE LEFT_JUSTIFY(STRING)

C  Thomas Matejka NOAA/NSSL 17 May 2000

      IMPLICIT NONE
      CHARACTER(LEN=*)::STRING
      INTEGER::I,L

      IF(STRING.NE.'')THEN
         L=LEN(STRING)
         DO I=1,L
            IF(STRING(I:I).NE.'')THEN
               EXIT
            ENDIF
         ENDDO
         IF(I.GT.1)THEN
            STRING(1:L-I+1)=STRING(I:L)
            STRING(L-I+2:L)=''
         ENDIF
      ENDIF

      END SUBROUTINE LEFT_JUSTIFY
