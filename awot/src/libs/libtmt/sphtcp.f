      REAL FUNCTION SPHTCP(TD,P)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL SHUM
      REAL TD,P,Q

      Q=SHUM(TD,P)
      SPHTCP=(1.-Q)*CP_DRY+Q*CP_VAP
      RETURN
      END
