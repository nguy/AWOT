      REAL FUNCTION AHUMS(T)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL VAPPS
      REAL T

      AHUMS=MW_H2O*VAPPS(T)/UNIV_GAS/T
      RETURN
      END
