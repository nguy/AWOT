      SUBROUTINE SF_I(I,F)

C  Thomas Matejka NOAA/NSSL 21 October 1998

      IMPLICIT NONE
      CHARACTER(LEN=*)::F
      INTEGER::I,N

      IF(I.EQ.0)THEN
         N=1
      ELSE
         N=IFIX(ALOG10(FLOAT(IABS(I))))+1
         IF(I.LT.0)THEN
            N=N+1
         ENDIF
      ENDIF

      F=''
      CALL APPEND_STRING(0,'I',F)
      CALL APPEND_INTEGER(0,N,F)

      END SUBROUTINE SF_I
