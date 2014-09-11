      REAL FUNCTION SHUMSI(T,P)

C  Thomas Matejka NOAA/NSSL 10 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL VAPPSI
      REAL T,P,ESI

      IF(T.LE.T_FREEZE)THEN
         ESI=VAPPSI(T)
         IF(ESI.GE.P)THEN
            WRITE(TMTLIB_MESSAGE_UNIT,*)'SHUMSI:  ESI IS GREATER THAN ',
     $      'OR EQUAL TO P.'
            STOP
         ENDIF
         SHUMSI=ESI*MW_H2O/((P-ESI)*MW_DRY+ESI*MW_H2O)
      ELSE
         SHUMSI=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
