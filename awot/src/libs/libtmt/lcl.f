      SUBROUTINE LCL(T0,TD0,P0,TLCL,PLCL)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPP,AMIX,AMIXS
      LOGICAL AGAIN
      INTEGER NITER
      REAL T0,TD0,P0,TLCL,PLCL,E0,ROVCP,W0,PUPPER,PLOWER,WSLCL

      E0=VAPP(TD0)
      IF(E0.GE.P0)THEN
         WRITE(TMTLIB_MESSAGE_UNIT,*)'LCL:  E0 IS GREATER THAN OR ',
     $   'EQUAL TO P0.'
         STOP
      ENDIF
      ROVCP=UNIV_GAS*P0/(CP_DRY*(P0-E0)*MW_DRY+CP_VAP*E0*MW_H2O)
      W0=AMIX(TD0,P0)
      IF(W0.LT.AMIXS(T0,P0)-WTOL)THEN
         PUPPER=0.
         PLOWER=P0
         AGAIN=.TRUE.
         NITER=0
         DOWHILE(AGAIN)
            NITER=NITER+1
            IF(NITER.GT.MAX_ITERS)THEN
               WRITE(TMTLIB_MESSAGE_UNIT,*)'LCL:  EXCEEDED MAXIMUM ',
     $         'ITERATIONS.'
               STOP
            ENDIF
            PLCL=(PUPPER+PLOWER)/2.
            TLCL=T0*(PLCL/P0)**ROVCP
            WSLCL=AMIXS(TLCL,PLCL)
            IF(WSLCL.GT.W0+WTOL)THEN
               PLOWER=PLCL
            ELSEIF(WSLCL.LT.W0-WTOL)THEN
               PUPPER=PLCL
            ELSE
               AGAIN=.FALSE.
            ENDIF
         ENDDO
      ELSE
         TLCL=T0
         PLCL=P0
      ENDIF
      RETURN
      END
