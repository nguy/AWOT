      REAL FUNCTION DPDZHS(T,TD,P)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL DEN
      REAL T,TD,P

      DPDZHS=-DEN(T,TD,P)*E_GRAV
      RETURN
      END
