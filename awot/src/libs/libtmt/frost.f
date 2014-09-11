      REAL FUNCTION FROST(TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL AMIX,AMIXSI
      LOGICAL AGAIN
      INTEGER NITER
      REAL TD,P,W,TFP1,TFP2,WSI

      W=AMIX(TD,P)
      IF(AMIXSI(T_FREEZE,P).GE.W)THEN
         IF(TD.GT.0.)THEN
            TFP1=0.
            TFP2=T_FREEZE
            AGAIN=.TRUE.
            NITER=0
            DOWHILE(AGAIN)
               NITER=NITER+1
               IF(NITER.GT.MAX_ITERS)THEN
                  WRITE(TMTLIB_MESSAGE_UNIT,*)'FROST:  EXCEEDED ',
     $            'MAXIMUM ITERATIONS.'
                  STOP
               ENDIF
               FROST=(TFP1+TFP2)/2.
               WSI=AMIXSI(FROST,P)
               IF(WSI.LT.W-WTOL)THEN
                  TFP1=FROST
               ELSEIF(WSI.GT.W+WTOL)THEN
                  TFP2=FROST
               ELSE
                  AGAIN=.FALSE.
               ENDIF
            ENDDO
         ELSE
            FROST=0.
         ENDIF
      ELSE
         FROST=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
