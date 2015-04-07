      REAL FUNCTION REFLECTIVITY_FACTOR(POWER_T,GAIN,BW_AZ,BW_EL,
     $PULSE_DUR,WAVELENGTH,K2,RANGE,POWER_R)

C  Thomas Matejka NOAA/NSSL 22 March 1993

C  This function returns the equivalent reflectivity factor of a volume
C  target.

C  POWER_T is the transmitted power.

C  GAIN is the antenna gain.

C  BW_AZ is the beam width in the azimuthal direction (deg).

C  BW_EL is the beam width in the elevation direction (deg).

C  PULSE_DUR is the pulse duration.

C  WAVELENGTH is the wavelength.

C  K2 is the term |K|**2 involving the index of refraction of the
C  targets.

C  RANGE is the distance to the target.

C  POWER_R is the echo power received at the radar.

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL POWER_T,GAIN,BW_AZ,BW_EL,PULSE_DUR,WAVELENGTH,K2,RANGE,
     $POWER_R,BW_AZ_RAD,BW_EL_RAD

      BW_AZ_RAD=BW_AZ*RADDEG
      BW_EL_RAD=BW_EL*RADDEG

      REFLECTIVITY_FACTOR=1024.*ALOG(2.)*WAVELENGTH**2*RANGE**2*POWER_R/
     $(C_LIGHT*PI**3*POWER_T*GAIN**2*BW_AZ_RAD*BW_EL_RAD*PULSE_DUR*K2)
      RETURN
      END
