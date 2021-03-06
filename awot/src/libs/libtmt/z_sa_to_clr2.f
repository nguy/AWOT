      REAL FUNCTION Z_SA_TO_CLR2(Z1,TV1,Z2,TV2,Z3,DEN3,Z_SA)

C  Thomas Matejka NOAA/NSSL 18 November 1997

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL,EXTERNAL::P_SA,ZCLR,TVCLR
      REAL::Z1,TV1,Z2,TV2,Z3,P3,Z_SA,P,TV3,DEN3

      P=P_SA(Z_SA)
      TV3=TVCLR(Z1,TV1,Z2,TV2,Z3)
      P3=DEN3*UNIV_GAS*TV3/MW_DRY
      Z_SA_TO_CLR2=ZCLR(Z1,TV1,Z2,TV2,Z3,P3,P)
      END FUNCTION Z_SA_TO_CLR2
