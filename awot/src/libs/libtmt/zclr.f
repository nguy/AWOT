      REAL FUNCTION ZCLR(Z1,TV1,Z2,TV2,Z3,P3,P)

C  Thomas Matejka NOAA/NSSL 18 November 1997

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL::Z1,TV1,Z2,TV2,Z3,P3,GAMMA,P

      IF(Z2.NE.Z1)THEN
         IF(TV2.EQ.TV1)THEN
            ZCLR=Z3-((UNIV_GAS*TV1)/(MW_DRY*E_GRAV))*ALOG(P/P3)
         ELSE
            GAMMA=(TV2-TV1)/(Z2-Z1)
            ZCLR=Z1+((TV1+GAMMA*(Z3-Z1))*(P/P3)**
     $      (-((UNIV_GAS*GAMMA)/(MW_DRY*E_GRAV))-TV1)/GAMMA)
         ENDIF
      ELSE
         WRITE(TMTLIB_MESSAGE_UNIT,*)'PCLR:  Z1 AND Z2 MUST NOT BE ',
     $   'EQUAL.'
         STOP
      ENDIF
      END FUNCTION ZCLR
