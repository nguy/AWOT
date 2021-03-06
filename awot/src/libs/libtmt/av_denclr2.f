      REAL FUNCTION AV_DENCLR2(Z1,TV1,Z2,TV2,Z3,P3,Z1_LAYER,Z2_LAYER)

C  Thomas Matejka NOAA/NSSL 20 July 1994

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL DENCLR2,TVCLR
      REAL Z1,TV1,Z2,TV2,Z3,P3,GAMMA,Z1_LAYER,Z2_LAYER,ALPHA,BETA,TV3,
     $DEN3

      IF(Z1_LAYER.NE.Z2_LAYER)THEN
         IF(Z2.NE.Z1)THEN
            TV3=TVCLR(Z1,TV1,Z2,TV2,Z3)
            DEN3=P3*MW_DRY/UNIV_GAS/TV3            
            IF(TV2.EQ.TV1)THEN
               BETA=-E_GRAV*MW_DRY/(UNIV_GAS*TV1)
               AV_DENCLR2=DEN3*(EXP(BETA*(Z2_LAYER-Z3))-
     $         EXP(BETA*(Z1_LAYER-Z3)))/((Z2_LAYER-Z1_LAYER)*BETA)
            ELSE
               GAMMA=(TV2-TV1)/(Z2-Z1)
               ALPHA=-1.-E_GRAV*MW_DRY/(UNIV_GAS*GAMMA)
               AV_DENCLR2=DEN3*(DBLE(TV1+GAMMA*(Z2_LAYER-Z1))**
     $         (ALPHA+1.)-DBLE(TV1+GAMMA*(Z1_LAYER-Z1))**(ALPHA+1.))/
     $         ((Z2_LAYER-Z1_LAYER)*DBLE(TV1+GAMMA*(Z3-Z1))**ALPHA*
     $         (ALPHA+1.)*GAMMA)
            ENDIF
         ELSE
            WRITE(TMTLIB_MESSAGE_UNIT,*)'AV_DENCLR2:  Z1 AND Z2 MUST ',
     $      'NOT BE EQUAL.'
            STOP
         ENDIF
      ELSE
         AV_DENCLR2=DENCLR2(Z1,TV1,Z2,TV2,Z3,P3,Z1_LAYER)
      ENDIF
      RETURN
      END
