      REAL FUNCTION SPGASC(TD,P)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL SHUM
      REAL TD,P,Q

      Q=SHUM(TD,P)
      SPGASC=UNIV_GAS*((1.-Q)/MW_DRY+Q/MW_H2O)
      RETURN
      END
