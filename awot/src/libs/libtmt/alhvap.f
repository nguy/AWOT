      REAL FUNCTION ALHVAP(T)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL T

      ALHVAP=LH_VAP_3PT+(CP_VAP-C_WAT)*(T-T_3PT)
      RETURN
      END
