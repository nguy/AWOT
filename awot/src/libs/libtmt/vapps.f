      REAL FUNCTION VAPPS(T)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL T,A,B

      IF(T.NE.0.)THEN
         A=(C_WAT-CP_VAP)*MW_H2O/UNIV_GAS
         B=LH_VAP_3PT*MW_H2O/UNIV_GAS/T_3PT
         VAPPS=P_3PT*(T_3PT/T)**A*EXP((A+B)*(1.-T_3PT/T))
      ELSE
         VAPPS=0.
      ENDIF
      RETURN
      END
