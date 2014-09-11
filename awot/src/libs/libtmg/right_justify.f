      SUBROUTINE RIGHT_JUSTIFY(STRING)

C  Thomas Matejka NOAA/NSSL 17 May 2000

      IMPLICIT NONE
      CHARACTER(LEN=*)::STRING
      INTEGER::I,L

      IF(STRING.NE.'')THEN
         L=LEN(STRING)
         DO I=L,1,-1
            IF(STRING(I:I).NE.'')THEN
               EXIT
            ENDIF
         ENDDO
         IF(I.LT.L)THEN
            STRING(L-I+1:L)=STRING(1:I)
            STRING(1:L-I)=''
         ENDIF
      ENDIF

      END SUBROUTINE RIGHT_JUSTIFY
