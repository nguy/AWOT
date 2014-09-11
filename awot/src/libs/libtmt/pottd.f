      REAL FUNCTION POTTD(T,TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPP
      REAL T,TD,P,E

      E=VAPP(TD)
      IF(E.GE.P)THEN
         WRITE(TMTLIB_MESSAGE_UNIT,*)'POTTD:  E IS GREATER THAN OR ',
     $   'EQUAL TO P.'
         STOP
      ENDIF
      POTTD=T*(P_REF/(P-E))**(UNIV_GAS/MW_DRY/CP_DRY)
      RETURN
      END
