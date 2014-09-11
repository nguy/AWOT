      REAL FUNCTION DPRH(T,P,U)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPPS,DPMIX
      REAL T,P,U,E,W

      E=U*VAPPS(T)
      IF(E.GE.P)THEN
         WRITE(TMTLIB_MESSAGE_UNIT,*)'DPRH:  E IS GREATER THAN OR ',
     $   'EQUAL TO P.'
         STOP
      ENDIF
      W=E*MW_H2O/(P-E)/MW_DRY
      DPRH=DPMIX(P,W)
      RETURN
      END
