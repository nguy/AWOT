      SUBROUTINE LDL(T0,TD0,P0,TLDL,PLDL)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPP,AMIX,AMIXSI
      LOGICAL AGAIN
      INTEGER NITER
      REAL T0,TD0,P0,TLDL,PLDL,E0,ROVCP,W0,PUPPER,PLOWER,WSILDL,
     $WSI_FREEZE,P_FREEZE

      E0=VAPP(TD0)
      IF(E0.GE.P0)THEN
         WRITE(TMTLIB_MESSAGE_UNIT,*)'LDL:  E0 IS GREATER THAN OR ',
     $   'EQUAL TO P0.'
         STOP
      ENDIF
      ROVCP=UNIV_GAS*P0/(CP_DRY*(P0-E0)*MW_DRY+CP_VAP*E0*MW_H2O)
      W0=AMIX(TD0,P0)
      P_FREEZE=P0*(T_FREEZE/T0)**(1./ROVCP)
      IF(T0.GT.T_FREEZE)THEN
         WSI_FREEZE=AMIXSI(T_FREEZE,P_FREEZE)
         IF(W0.LT.WSI_FREEZE-WTOL)THEN
            PUPPER=0.
            PLOWER=P_FREEZE
            AGAIN=.TRUE.
            NITER=0
            DOWHILE(AGAIN)
               NITER=NITER+1
               IF(NITER.GT.MAX_ITERS)THEN
                  WRITE(TMTLIB_MESSAGE_UNIT,*)'LDL:  EXCEEDED MAXIMUM ',
     $            'ITERATIONS.'
                  STOP
               ENDIF
               PLDL=(PUPPER+PLOWER)/2.
               TLDL=T0*(PLDL/P0)**ROVCP
               WSILDL=AMIXSI(TLDL,PLDL)
               IF(WSILDL.GT.W0+WTOL)THEN
                  PLOWER=PLDL
               ELSEIF(WSILDL.LT.W0-WTOL)THEN
                  PUPPER=PLDL
               ELSE
                  AGAIN=.FALSE.
               ENDIF
            ENDDO
         ELSEIF(W0.GT.WSI_FREEZE)THEN
            TLDL=TMTLIB_BADFLAG
            PLDL=TMTLIB_BADFLAG
         ELSE
            TLDL=T_FREEZE
            PLDL=P_FREEZE
         ENDIF
      ELSE
         IF(W0.LT.AMIXSI(T0,P0)-WTOL)THEN
            PUPPER=0.
            PLOWER=P_FREEZE
            AGAIN=.TRUE.
            NITER=0
            DOWHILE(AGAIN)
               NITER=NITER+1
               IF(NITER.GT.MAX_ITERS)THEN
                  WRITE(TMTLIB_MESSAGE_UNIT,*)'LDL:  EXCEEDED MAXIMUM ',
     $            'ITERATIONS.'
                  STOP
               ENDIF
               PLDL=(PUPPER+PLOWER)/2.
               TLDL=T0*(PLDL/P0)**ROVCP
               WSILDL=AMIXSI(TLDL,PLDL)
               IF(WSILDL.GT.W0+WTOL)THEN
                  PLOWER=PLDL
               ELSEIF(WSILDL.LT.W0-WTOL)THEN
                  PUPPER=PLDL
               ELSE
                  AGAIN=.FALSE.
               ENDIF
            ENDDO
         ELSE
            TLDL=T0
            PLDL=P0
         ENDIF
      ENDIF
      RETURN
      END
