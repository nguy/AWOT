      REAL FUNCTION DPFP(P,TF)

C  Thomas Matejka NOAA/NSSL 9 March 1993

      IMPLICIT NONE
      INCLUDE 'tmtlib.inc'
      REAL AMIXSI,DPMIX
      REAL P,TF,W

      W=AMIXSI(TF,P)
      IF(W.NE.TMTLIB_BADFLAG)THEN
         DPFP=DPMIX(P,W)
      ELSE
         DPFP=TMTLIB_BADFLAG
      ENDIF
      RETURN
      END
