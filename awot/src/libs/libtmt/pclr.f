      REAL FUNCTION PCLR(Z1,TV1,Z2,TV2,Z3,P3,Z)

C  Thomas Matejka NOAA/NSSL 9 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL Z1,TV1,Z2,TV2,Z3,P3,Z,GAMMA

      IF(Z2.NE.Z1)THEN
         IF(TV2.EQ.TV1)THEN
            PCLR=P3*EXP(-E_GRAV*(Z-Z3)*MW_DRY/UNIV_GAS/TV1)
         ELSE
            GAMMA=(TV2-TV1)/(Z2-Z1)
            PCLR=P3*((TV1+GAMMA*(Z-Z1))/(TV1+GAMMA*(Z3-Z1)))**
     $      (-(E_GRAV*MW_DRY/UNIV_GAS/GAMMA))
         ENDIF
      ELSE
         WRITE(TMTLIB_MESSAGE_UNIT,*)'PCLR:  Z1 AND Z2 MUST NOT BE ',
     $   'EQUAL.'
         STOP
      ENDIF
      RETURN
      END
