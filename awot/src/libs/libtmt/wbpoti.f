      REAL FUNCTION WBPOTI(T,TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      LOGICAL AGAIN
      INTEGER NINC
      REAL T,TD,P,TLDL,PLDL,P2,T2,DELP,P1,T1,PAVE,DTDP,DQSDT,DQSDP,TAVE

      CALL LDL(T,TD,P,TLDL,PLDL)
      IF(PLDL.NE.TMTLIB_BADFLAG)THEN
         P2=PLDL
         T2=TLDL
         DELP=SIGN(DELTA_P,P_REF-PLDL)
         AGAIN=.TRUE.
         NINC=0
         DOWHILE(AGAIN)
            NINC=NINC+1
            IF(NINC.GT.MAX_INCS)THEN
               WRITE(TMTLIB_MESSAGE_UNIT,*)'WBPOTI:  EXCEEDED MAXIMUM ',
     $         'INCREMENTS.'
               STOP
            ENDIF
            P1=P2
            T1=T2
            IF(ABS(P1-P_REF).LE.ABS(DELP))THEN
               DELP=P_REF-P1
               AGAIN=.FALSE.
            ENDIF
            P2=P1+DELP
            PAVE=(P1+P2)/2.
            IF(T1.LT.T_FREEZE)THEN
               CALL PSEUDI(T1,PAVE,DTDP,DQSDT,DQSDP)
            ELSE
               CALL PSEUDO(T1,PAVE,DTDP,DQSDT,DQSDP)
            ENDIF
            T2=T1+DTDP*DELP
            TAVE=(T1+T2)/2.
            IF(TAVE.LT.T_FREEZE)THEN
               CALL PSEUDI(TAVE,PAVE,DTDP,DQSDT,DQSDP)
            ELSE
               CALL PSEUDO(TAVE,PAVE,DTDP,DQSDT,DQSDP)
            ENDIF
            T2=T1+DTDP*DELP
            IF(T1.LT.T_FREEZE.AND.T2.GT.T_FREEZE)THEN
               AGAIN=.TRUE.
               T2=T_FREEZE
               TAVE=(T1+T2)/2.
               CALL PSEUDI(TAVE,P1,DTDP,DQSDT,DQSDP)
               P2=P1+(T2-T1)/DTDP
               PAVE=(P1+P2)/2.
               CALL PSEUDI(TAVE,PAVE,DTDP,DQSDT,DQSDP)
               P2=P1+(T2-T1)/DTDP
            ELSEIF(T1.GT.T_FREEZE.AND.T2.LT.T_FREEZE)THEN
               AGAIN=.TRUE.
               T2=T_FREEZE
               TAVE=(T1+T2)/2.
               CALL PSEUDO(TAVE,P1,DTDP,DQSDT,DQSDP)
               P2=P1+(T2-T1)/DTDP
               PAVE=(P1+P2)/2.
               CALL PSEUDO(TAVE,P1,DTDP,DQSDT,DQSDP)
               P2=P1+(T2-T1)/DTDP
            ENDIF
         ENDDO
         WBPOTI=T2
      ELSE
         WBPOTI=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
