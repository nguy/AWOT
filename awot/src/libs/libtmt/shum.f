      REAL FUNCTION SHUM(TD,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPP
      REAL TD,P,E

      E=VAPP(TD)
      IF(E.GE.P)THEN
         WRITE(TMTLIB_MESSAGE_UNIT,*)'SHUM:  E IS GREATER THAN OR ',
     $   'EQUAL TO P.'
         STOP
      ENDIF
      SHUM=E*MW_H2O/((P-E)*MW_DRY+E*MW_H2O)
      RETURN
      END
