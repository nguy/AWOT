      REAL FUNCTION AHUMSI(T)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPPSI
      REAL T

      IF(T.LE.T_FREEZE)THEN
         AHUMSI=MW_H2O*VAPPSI(T)/UNIV_GAS/T
      ELSE
         AHUMSI=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
