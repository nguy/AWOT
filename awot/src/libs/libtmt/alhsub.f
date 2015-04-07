      REAL FUNCTION ALHSUB(T)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL T

      IF(T.LE.T_FREEZE)THEN
         ALHSUB=LH_SUB_3PT+(CP_VAP-C_ICE)*(T-T_3PT)
      ELSE
         ALHSUB=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
