      REAL FUNCTION DPRHI(T,P,UI)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPPSI,DPMIX
      REAL T,P,UI,ESI,E,W

      ESI=VAPPSI(T)
      IF(ESI.NE.TMTLIB_BADFLAG)THEN
         E=UI*ESI
         IF(E.GE.P)THEN
            WRITE(TMTLIB_MESSAGE_UNIT,*)'DPRHI:  E IS GREATER THAN OR ',
     $      'EQUAL TO P.'
            STOP
         ENDIF
         W=E*MW_H2O/(P-E)/MW_DRY
         DPRHI=DPMIX(P,W)
      ELSE
         DPRHI=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
