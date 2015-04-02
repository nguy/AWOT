      REAL FUNCTION POTT(T,TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPP
      REAL T,TD,P,E

      E=VAPP(TD)
      IF(E.GE.P)THEN
         WRITE(TMTLIB_MESSAGE_UNIT,*)'POTT:  E IS GREATER THAN OR ',
     $   'EQUAL TO P.'
         STOP
      ENDIF
      POTT=T*(P_REF/P)**(UNIV_GAS*P/(CP_DRY*(P-E)*MW_DRY+
     $CP_VAP*E*MW_H2O))
      RETURN
      END
