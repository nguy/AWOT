      REAL FUNCTION Z_SA_TO_CLR(Z1,TV1,Z2,TV2,Z3,P3,Z_SA)

C  Thomas Matejka NOAA/NSSL 18 November 1997

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL,EXTERNAL::P_SA,ZCLR
      REAL::Z1,TV1,Z2,TV2,Z3,P3,Z_SA,P

      P=P_SA(Z_SA)
      Z_SA_TO_CLR=ZCLR(Z1,TV1,Z2,TV2,Z3,P3,P)
      END FUNCTION Z_SA_TO_CLR
