      SUBROUTINE SF_R_TRUNC(A,D,F)

C  Thomas Matejka NOAA/NSSL 21 October 1998

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER(LEN=MAX_STRING)::STRING
      CHARACTER(LEN=*)::F
      INTEGER::N,D,I
      REAL::A,B

      IF(D.LT.0)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'SF_R_TRUNC:  D MUST BE AT LEAST ',
     $   '0.'
         STOP
      ENDIF

      IF(A.EQ.0.)THEN
         F='I1'
         RETURN
      ENDIF

      B=ABS(A)
      IF(B.GE.1.)THEN
         N=IFIX(ALOG10(B))+2+D
      ELSE
         N=1+D
      ENDIF
      IF(A.LT.0.)THEN
         N=N+1
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
            CALL APPEND_INTEGER(0,I-1,F)
            EXIT
         ELSEIF(STRING(I:I).NE.'0')THEN
            F=''
            CALL APPEND_STRING(0,'F',F)
            CALL APPEND_INTEGER(0,I,F)
            CALL APPEND_STRING(0,'.',F)
            CALL APPEND_INTEGER(0,D-N+I,F)
            EXIT
         ENDIF
      ENDDO

      END SUBROUTINE SF_R_TRUNC
