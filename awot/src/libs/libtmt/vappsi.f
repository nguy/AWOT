      REAL FUNCTION VAPPSI(T)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL T,A,B

      IF(T.LE.T_FREEZE)THEN
         IF(T.NE.0.)THEN
            A=(C_ICE-CP_VAP)*MW_H2O/UNIV_GAS
            B=LH_SUB_3PT*MW_H2O/UNIV_GAS/T_3PT
            VAPPSI=P_3PT*(T_3PT/T)**A*EXP((A+B)*(1.-T_3PT/T))
         ELSE
            VAPPSI=0.
         ENDIF
      ELSE
         VAPPSI=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
