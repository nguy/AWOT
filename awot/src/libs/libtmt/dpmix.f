      REAL FUNCTION DPMIX(P,W)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'tmtlib.inc'
      REAL AMIXS
      LOGICAL AGAIN
      INTEGER NITER
      REAL P,W,TD1,TD2,WS

      IF(W.GT.0.)THEN
         TD1=0.
         TD2=TD_MAX
         AGAIN=.TRUE.
         NITER=0
         DOWHILE(AGAIN)
            NITER=NITER+1
            IF(NITER.GT.MAX_ITERS)THEN
               WRITE(TMTLIB_MESSAGE_UNIT,*)'DPMIX:  EXCEEDED MAXIMUM ',
     $         'ITERATIONS.'
               STOP
            ENDIF
            DPMIX=(TD1+TD2)/2.
            WS=AMIXS(DPMIX,P)
            IF(WS.LT.W-WTOL)THEN
               TD1=DPMIX
            ELSEIF(WS.GT.W+WTOL)THEN
               TD2=DPMIX
            ELSE
               AGAIN=.FALSE.
            ENDIF
         ENDDO
      ELSE
         DPMIX=0.
      ENDIF
      RETURN
      END
