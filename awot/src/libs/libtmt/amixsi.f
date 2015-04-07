      REAL FUNCTION AMIXSI(T,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPPSI
      REAL T,P,ESI

      IF(T.LE.T_FREEZE)THEN
         ESI=VAPPSI(T)
         IF(ESI.GE.P)THEN
            WRITE(TMTLIB_MESSAGE_UNIT,*)'AMIXSI:  ESI IS GREATER THAN ',
     $      'OR EQUAL TO P.'
            STOP
         ENDIF
         AMIXSI=ESI*MW_H2O/(P-ESI)/MW_DRY
      ELSE
         AMIXSI=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
