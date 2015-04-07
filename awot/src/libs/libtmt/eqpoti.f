      REAL FUNCTION EQPOTI(T,TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL POTTD,ALHSUB,AMIXSI
      LOGICAL AGAIN
      INTEGER NINC
      REAL T,TD,P,SUM,TLDL,PLDL,THETAD,ALSLDL,WSILDL,T2,P2,WSI2,T1,P1,
     $WSI1,TAVE,PAVE,DTDPI,DQSIDTI,DQSIDPI,WSIAVE,DELT

      DELT=-DELTA_T
      SUM=0.
      CALL LDL(T,TD,P,TLDL,PLDL)
      IF(TLDL.NE.TMTLIB_BADFLAG)THEN
         THETAD=POTTD(TLDL,TLDL,PLDL)
         ALSLDL=ALHSUB(TLDL)
         WSILDL=AMIXSI(TLDL,PLDL)
         T2=TLDL
         P2=PLDL
         WSI2=WSILDL
         AGAIN=.TRUE.
         NINC=0
         DOWHILE(AGAIN)
            NINC=NINC+1
            IF(NINC.GT.MAX_INCS)THEN
               WRITE(TMTLIB_MESSAGE_UNIT,*)'EQPOTI:  EXCEEDED MAXIMUM ',
     $         'INCREMENTS.'
               STOP
            ENDIF
            T1=T2
            P1=P2
            WSI1=WSI2
            DOWHILE(T1+DELT.LE.0.)
               DELT=DELT/2.
            ENDDO
            T2=T1+DELT
            TAVE=(T1+T2)/2.
            CALL PSEUDI(TAVE,P1,DTDPI,DQSIDTI,DQSIDPI)
            P2=P1+DELT/DTDPI
            PAVE=(P1+P2)/2.
            CALL PSEUDI(TAVE,PAVE,DTDPI,DQSIDTI,DQSIDPI)
            P2=P1+DELT/DTDPI
            WSI2=AMIXSI(T2,P2)
            WSIAVE=(WSI1+WSI2)/2.
            SUM=SUM+WSIAVE*C_ICE*ALOG(T2/T1)/CP_DRY
            IF(WSI2.LE.WTOL)THEN
               AGAIN=.FALSE.
            ENDIF
         ENDDO
         EQPOTI=THETAD*EXP(ALSLDL*WSILDL/TLDL/CP_DRY-SUM)
      ELSE
         EQPOTI=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
