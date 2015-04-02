      REAL FUNCTION PCLR2(Z1,TV1,Z2,TV2,Z3,DEN3,Z)

C  Thomas Matejka NOAA/NSSL 10 May 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL TVCLR
      REAL Z1,TV1,Z2,TV2,Z3,DEN3,TV3,P3,Z,GAMMA

      IF(Z2.NE.Z1)THEN
         TV3=TVCLR(Z1,TV1,Z2,TV2,Z3)
         P3=DEN3*UNIV_GAS*TV3/MW_DRY
         IF(TV2.EQ.TV1)THEN
            PCLR2=P3*EXP(-E_GRAV*(Z-Z3)*MW_DRY/UNIV_GAS/TV1)
         ELSE
            GAMMA=(TV2-TV1)/(Z2-Z1)
            PCLR2=P3*((TV1+GAMMA*(Z-Z1))/(TV1+GAMMA*(Z3-Z1)))**
     $      (-(E_GRAV*MW_DRY/UNIV_GAS/GAMMA))
         ENDIF
      ELSE
         WRITE(TMTLIB_MESSAGE_UNIT,*)'PCLR2:  Z1 AND Z2 MUST NOT BE ',
     $   'EQUAL.'
         STOP
      ENDIF
      RETURN
      END
