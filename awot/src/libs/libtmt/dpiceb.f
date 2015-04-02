      REAL FUNCTION DPICEB(T,P,TI)

C  Thomas Matejka NOAA/NSSL 8 July 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      INCLUDE 'tmtlib.inc'
      REAL ALHSUB,AMIXSI,DPMIX
      REAL T,P,TI,ALS,W

      IF(TI.GT.0.)THEN
         ALS=ALHSUB(TI)
         IF(ALS.NE.TMTLIB_BADFLAG)THEN
            W=(AMIXSI(TI,P)*ALS-(T-TI)*CP_DRY)/((T-TI)*CP_VAP+ALS)
            DPICEB=DPMIX(P,W)
         ELSE
            DPICEB=TMTLIB_BADFLAG
         ENDIF
      ELSE
         DPICEB=0.
      ENDIF
      RETURN
      END
