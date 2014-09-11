      REAL FUNCTION DPSHUM(P,Q)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'tmtlib.inc'
      REAL DPMIX
      REAL P,Q,W

      IF(Q.GE.1.)THEN
         WRITE(TMTLIB_MESSAGE_UNIT,*)'DPSHUM:  Q IS GREATER THAN OR ',
     $   'EQUAL TO 1..'
         STOP
      ENDIF
      W=Q/(1.-Q)
      DPSHUM=DPMIX(P,W)
      RETURN
      END
