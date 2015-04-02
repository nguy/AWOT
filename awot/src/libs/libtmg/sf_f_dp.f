      SUBROUTINE SF_R_DP(A,D,F)

C  Thomas Matejka NOAA/NSSL 17 May 2000

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      CHARACTER(LEN=*)::F
      INTEGER::N,D
      REAL::A,B

      IF(D.LT.0)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'SF_R_DP:  D MUST BE AT LEAST 0.'
         STOP
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

      END SUBROUTINE SF_R_DP
