      REAL FUNCTION WETBLB(T,TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL AMIX,AMIXS,ALHVAP
      LOGICAL AGAIN
      INTEGER NITER
      REAL T,TD,P,W0,TW1,TW2,W,TGUESS

      W0=AMIX(TD,P)
      TW1=0.
      TW2=T_MAX
      AGAIN=.TRUE.
      NITER=0
      DOWHILE(AGAIN)
         NITER=NITER+1
         IF(NITER.GT.MAX_ITERS)THEN
            WRITE(TMTLIB_MESSAGE_UNIT,*)'WETBLB:  EXCEEDED MAXIMUM ',
     $      'ITERATIONS.'
            STOP
         ENDIF
         WETBLB=(TW1+TW2)/2.
         W=AMIXS(WETBLB,P)
         TGUESS=WETBLB+(W-W0)*ALHVAP(WETBLB)/(CP_DRY+W0*CP_VAP)
         IF(TGUESS.LT.T-TTOL)THEN
            TW1=WETBLB
         ELSEIF(TGUESS.GT.T+TTOL)THEN
            TW2=WETBLB
         ELSE
            AGAIN=.FALSE.
         ENDIF
      ENDDO
      RETURN
      END
