      SUBROUTINE SF_R_FIX_TRUNC(A,N,D,F)

C  Thomas Matejka NOAA/NSSL 21 October 1998

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER(LEN=MAX_STRING)::STRING
      CHARACTER(LEN=*)::F
      INTEGER::N,D,I
      REAL::A

      IF(N.LT.1)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'SF_R_FIX_TRUNC:  N MUST BE AT ',
     $   'LEAST 1.'
         STOP
      ENDIF
      IF(D.LT.0)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'SF_R_FIX_TRUNC:  D MUST BE AT ',
     $   'LEAST 0.'
         STOP
      ENDIF

      IF(A.EQ.0.)THEN
         F=''
         CALL APPEND_STRING(0,'I',F)
         CALL APPEND_STRING(0,N,F)
         RETURN
      ENDIF

      F=''
      CALL APPEND_STRING(0,'F',F)
      CALL APPEND_INTEGER(0,N,F)
      CALL APPEND_STRING(0,'.',F)
      CALL APPEND_INTEGER(0,D,F)

      WRITE(STRING,"("//F//")")A
      DO I=N,1,-1
         IF(STRING(I:I).EQ.'.')THEN
            F=''
            CALL APPEND_STRING(0,'I',F)
            CALL APPEND_INTEGER(0,N,F)
            EXIT
         ELSEIF(STRING(I:I).NE.'0')THEN
            F=''
            CALL APPEND_STRING(0,'F',F)
            CALL APPEND_INTEGER(0,N,F)
            CALL APPEND_STRING(0,'.',F)
            CALL APPEND_INTEGER(0,D-N+I,F)
            EXIT
         ENDIF
      ENDDO

      END SUBROUTINE SF_R_FIX_TRUNC
