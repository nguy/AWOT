      REAL FUNCTION SPHTCV(TD,P)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL SHUM
      REAL TD,P,Q

      Q=SHUM(TD,P)
      SPHTCV=(1.-Q)*(CP_DRY-UNIV_GAS/MW_DRY)+Q*(CP_VAP-UNIV_GAS/MW_H2O)
      RETURN
      END
