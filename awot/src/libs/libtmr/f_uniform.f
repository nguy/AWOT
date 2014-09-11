      SUBROUTINE F_UNIFORM(Z_TOP_RAIN,Z_BOT_SNOW,F_RAIN,F_SNOW,Z,F,SD_F)

C  Thomas Matejka NOAA/NSSL 15 February 1995

C  This subroutine calculates the hydrometeor terminal fall speed and
C  its standard deviation from specified values of rain and snow
C  terminal fall speeds and altitude.

C  Terminal fall speeds for rain and snow are linearly interpolated in a
C  specified melting layer.

C  Formulas for the standard deviations of the terminal fall speeds are
C  purely empirical and attempt to take into account additional errors
C  in measured radar reflectivity factor.  The standard deviation is
C  boosted further in the melting layer.

C  Input:

C  Z_TOP_RAIN is a real variable that specifies the altitude (m) of the
C  top of the layer in which the hydrometeors are assumed to be liquid.

C  Z_BOT_SNOW is a real variable that specifies the altitude (m) of the
C  bottom of the layer in which the hydrometeors are assumed to be ice.
C  Z_BOT_SNOW must be greater than or equal to Z_TOP_RAIN.

C  F_RAIN is a real variable that specifies the terminal fall speed (m
C  s**-1) of liquid hydrometeors.

C  F_SNOW is a real variable that specifies the terminal fall speed (m
C  s**-1) of solid hydrometeors.

C  Z is a real variable that specifies the altitude (m) at which an
C  echo-power-weighted hydrometeor terminal fall speed estimate is
C  desired.

C  Output:

C  F is a real variable that returns an estimate of the hydrometeor
C  terminal fall speed (m s**-1).

C  SD_F is a real variable that returns the standard deviation of F (m
C  s**-1).

      IMPLICIT NONE
      REAL Z_TOP_RAIN,Z_BOT_SNOW,Z,F,SD_F,F_RAIN,F_SNOW,SD_F_RAIN,
     $SD_F_SNOW

C  Calculate the terminal fall speed and its standard deviation for
C  solid hydrometeors.
      IF(F_SNOW.LE.1.)THEN
         SD_F_SNOW=-0.5+1.*1.
      ELSE
         SD_F_SNOW=-0.5+1.*F_SNOW
      ENDIF

C  Calculate the terminal fall speed and its standard deviation for
C  liquid hydrometeors.
      SD_F_RAIN=0.6+0.2*F_RAIN

C  The hydrometeors are assumed to be ice.
      IF(Z.GE.Z_BOT_SNOW)THEN
         F=F_SNOW
         SD_F=SD_F_SNOW

C  The hydrometeors are assumed to be liquid.
      ELSEIF(Z.LE.Z_TOP_RAIN)THEN
         F=F_RAIN
         SD_F=SD_F_RAIN

C  The hydrometeors are assumed to be a mixture of ice and water.
C  Interpolate the values for ice and water and boost the standard
C  deviation.
      ELSE
         F=(F_RAIN*(Z_BOT_SNOW-Z)+F_SNOW*(Z-Z_TOP_RAIN))/
     $   (Z_BOT_SNOW-Z_TOP_RAIN)
         SD_F=(SD_F_RAIN*(Z_BOT_SNOW-Z)+SD_F_SNOW*(Z-Z_TOP_RAIN))/
     $   (Z_BOT_SNOW-Z_TOP_RAIN)
         SD_F=SD_F+1.0
      ENDIF

C  Done.
      RETURN
      END
