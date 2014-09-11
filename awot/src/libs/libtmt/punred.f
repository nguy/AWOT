      REAL FUNCTION PUNRED(T,TD,P,ZSTN,ZREF)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'tmtlib.inc'
      REAL PREDUC
      LOGICAL AGAIN
      INTEGER NITER
      REAL T,TD,P,ZSTN,ZREF,PMAX,PMIN,PTRIAL

      IF(ZSTN.EQ.ZREF)THEN
         PUNRED=P
      ELSE
         IF(ZSTN.GT.ZREF)THEN
            PMAX=P
            PMIN=0.
         ELSE
            PMAX=P_MAX
            PMIN=P
         ENDIF
         AGAIN=.TRUE.
         NITER=0
         DOWHILE(AGAIN)
            NITER=NITER+1
            IF(NITER.GT.MAX_ITERS)THEN
               WRITE(TMTLIB_MESSAGE_UNIT,*)'PUNRED:  EXCEEDED MAXIMUM ',
     $         'ITERATIONS.'
               STOP
            ENDIF
            PUNRED=(PMAX+PMIN)/2.
            PTRIAL=PREDUC(T,TD,PUNRED,ZSTN,ZREF)
            IF(PTRIAL.GT.P+PTOL)THEN
               PMAX=PUNRED
            ELSEIF(PTRIAL.LT.P-PTOL)THEN
               PMIN=PUNRED
            ELSE
               AGAIN=.FALSE.
            ENDIF
         ENDDO
      ENDIF
      RETURN
      END
