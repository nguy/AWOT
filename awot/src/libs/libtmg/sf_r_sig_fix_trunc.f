      SUBROUTINE SF_R_SIG_FIX_TRUNC(A,N,S,F)

C  Thomas Matejka NOAA/NSSL 21 October 1998

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER(LEN=MAX_STRING)::STRING
      CHARACTER(LEN=*)::F
      INTEGER::S,N,D,Q,I
      REAL::A,C

      IF(N.LT.1)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'SF_R_SIG_FIX_TRUNC:  N MUST BE ',
     $   'AT LEAST 1.'
         STOP
      ENDIF
      IF(S.LT.1)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'SF_R_SIG_FIX_TRUNC:  S MUST BE ',
     $   'AT LEAST 1.'
         STOP
      ENDIF

      IF(A.EQ.0.)THEN
         F=''
         CALL APPEND_STRING(0,'I',F)
         CALL APPEND_INTEGER(0,N,F)
         RETURN
      ENDIF

      C=ALOG10(ABS(A))
      Q=IFIX(C)
      IF(C.LT.0..AND.
     $AMOD(C,1.).NE.0.)THEN
         Q=Q-1
      ENDIF
      IF(Q.LT.0)THEN
         D=S-Q-1
      ELSEIF(Q.GE.S-1)THEN
         D=0
      ELSE
         D=S-Q-1
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

      END SUBROUTINE SF_R_SIG_FIX_TRUNC
