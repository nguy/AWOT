      REAL FUNCTION AHUM(T,TD)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL VAPP
      REAL T,TD

      AHUM=MW_H2O*VAPP(TD)/UNIV_GAS/T
      RETURN
      END
