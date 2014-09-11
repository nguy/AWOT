      REAL FUNCTION RHUMI(T,TD)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPP,VAPPSI
      REAL T,TD

      IF(T.LE.T_FREEZE)THEN
         RHUMI=VAPP(TD)/VAPPSI(T)
      ELSE
         RHUMI=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
