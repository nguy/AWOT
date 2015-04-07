      REAL FUNCTION WBPOTT(T,TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      LOGICAL AGAIN
      INTEGER NINC
      REAL T,TD,P,TLCL,PLCL,P2,T2,DELP,P1,T1,PAVE,DTDPP,DQSDTP,
     $DQSDPP,TAVE

      CALL LCL(T,TD,P,TLCL,PLCL)
      P2=PLCL
      T2=TLCL
      DELP=SIGN(DELTA_P,P_REF-PLCL)
      AGAIN=.TRUE.
      NINC=0
      DOWHILE(AGAIN)
         NINC=NINC+1
         IF(NINC.GT.MAX_INCS)THEN
            WRITE(TMTLIB_MESSAGE_UNIT,*)'WBPOTT:  EXCEEDED MAXIMUM ',
     $      'INCREMENTS.'
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
         CALL PSEUDO(T1,PAVE,DTDPP,DQSDTP,DQSDPP)
         T2=T1+DTDPP*DELP
         TAVE=(T1+T2)/2.
         CALL PSEUDO(TAVE,PAVE,DTDPP,DQSDTP,DQSDPP)
         T2=T1+DTDPP*DELP
      ENDDO
      WBPOTT=T2
      RETURN
      END
