      REAL FUNCTION VIRT(T,TD,P)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL SHUM
      REAL T,TD,P,Q

      Q=SHUM(TD,P)
      VIRT=T*(1.+Q*(MW_DRY/MW_H2O-1.))
      RETURN
      END
