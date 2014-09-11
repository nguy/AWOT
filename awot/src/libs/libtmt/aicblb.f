      REAL FUNCTION AICBLB(T,TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL AMIX,AMIXSI,ALHSUB
      LOGICAL AGAIN
      INTEGER NITER
      REAL T,TD,P,W0,TGUESS,TI1,TI2,W

      W0=AMIX(TD,P)
      TGUESS=T_FREEZE+(AMIXSI(T_FREEZE,P)-W0)*ALHSUB(T_FREEZE)/
     $(CP_DRY+W0*CP_VAP)
      IF(TGUESS.GE.T)THEN
         IF(TD.GT.0.)THEN
            TI1=0.
            TI2=T_FREEZE
            AGAIN=.TRUE.
            NITER=0
            DOWHILE(AGAIN)
               NITER=NITER+1
               IF(NITER.GT.MAX_ITERS)THEN
                  WRITE(TMTLIB_MESSAGE_UNIT,*)'AICBLB:  EXCEEDED ',
     $            'MAXIMUM ITERATIONS.'
                  STOP
               ENDIF
               AICBLB=(TI1+TI2)/2.
               W=AMIXSI(AICBLB,P)
               TGUESS=AICBLB+(W-W0)*ALHSUB(AICBLB)/(CP_DRY+W0*CP_VAP)
               IF(TGUESS.LT.T-TTOL)THEN
                  TI1=AICBLB
               ELSEIF(TGUESS.GT.T+TTOL)THEN
                  TI2=AICBLB
               ELSE
                  AGAIN=.FALSE.
               ENDIF
            ENDDO
         ELSE
            AICBLB=0.
         ENDIF
      ELSE
         AICBLB=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
