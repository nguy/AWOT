      REAL FUNCTION DENCLR(Z1,TV1,Z2,TV2,Z3,DEN3,Z)

C  Thomas Matejka NOAA/NSSL 9 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL Z1,TV1,Z2,TV2,Z3,DEN3,Z,GAMMA,ALPHA,BETA

      IF(Z2.NE.Z1)THEN
         IF(TV2.EQ.TV1)THEN
            BETA=-E_GRAV*MW_DRY/(UNIV_GAS*TV1)
            DENCLR=DEN3*EXP(BETA*(Z-Z3))
         ELSE
            GAMMA=(TV2-TV1)/(Z2-Z1)
            ALPHA=-1.-E_GRAV*MW_DRY/(UNIV_GAS*GAMMA)
            DENCLR=DEN3*((TV1+GAMMA*(Z-Z1))/(TV1+GAMMA*(Z3-Z1)))**ALPHA
         ENDIF
      ELSE
         WRITE(TMTLIB_MESSAGE_UNIT,*)'DENCLR:  Z1 AND Z2 MUST NOT BE ',
     $   'EQUAL.'
         STOP
      ENDIF
      RETURN
      END
