      REAL FUNCTION RHUM(T,TD)

C  Thomas Matejka NOAA/NSSL 8 March 1993

      IMPLICIT NONE
      REAL VAPP,VAPPS
      REAL T,TD

      RHUM=VAPP(TD)/VAPPS(T)
      RETURN
      END
