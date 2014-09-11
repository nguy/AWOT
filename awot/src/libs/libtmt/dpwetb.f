      REAL FUNCTION DPWETB(T,P,TW)

C  Thomas Matejka NOAA/NSSL 7 July 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL ALHVAP,AMIXS,DPMIX
      REAL T,P,TW,ALV,W

      IF(TW.GT.0.)THEN
         ALV=ALHVAP(TW)
         W=(AMIXS(TW,P)*ALV-(T-TW)*CP_DRY)/((T-TW)*CP_VAP+ALV)
         DPWETB=DPMIX(P,W)
      ELSE
         DPWETB=0.
      ENDIF
      RETURN
      END
