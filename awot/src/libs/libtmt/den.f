      REAL FUNCTION DEN(T,TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPP
      REAL T,TD,P,E

      E=VAPP(TD)
      IF(E.GE.P)THEN
         WRITE(TMTLIB_MESSAGE_UNIT,*)'DEN:  E IS GREATER THAN OR ',
     $   'EQUAL TO P.'
         STOP
      ENDIF
      DEN=((P-E)*MW_DRY+E*MW_H2O)/UNIV_GAS/T
      RETURN
      END
