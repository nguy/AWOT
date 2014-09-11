      SUBROUTINE F_FROM_DBZ(Z1,TV1,Z2,TV2,Z3,RHO3,Z_TOP_RAIN,Z_BOT_SNOW,
     $Z,DBZE,F,SD_F)

C  Thomas Matejka NOAA/NSSL 15 February 1995

C  This subroutine estimates the echo-power-weighted hydrometeor
C  terminal fall speed and its standard deviation from the reflectivity
C  factor and altitude.

C  The relation for rain assumes an exponential distribution with the
C  Marshall-Palmer (1948) value for N0 = 0.08 cm**-4 and the Gunn and
C  Kinzer (1949) terminal fall speed data, which includes asymptotic
C  behavior at large diameters (see Atlas et al., 1973).

C  The relation for snow is based on data of Gunn and Marshall (1958) as
C  interpreted by Sekhon and Srivastava (1970) (see Atlas et al., 1973).

C  Terminal fall speeds for rain and snow are linearly interpolated in a
C  specified melting layer.

C  The effect of air density is based on the complete formula of Foote
C  and Du Toit (1969).

C  Formulae for the standard deviations of the terminal fall speeds are
C  purely empirical and attempt to take into account additional errors
C  in measured radar reflectivity factor.  The standard deviation is
C  boosted further in the melting layer.

C  Input:

C  Z1 is a real variable that specifies an altitude (m).

C  TV1 is a real variable that specifies the virtual temperature (K) at
C  Z1.

C  Z2 is a real variable that specifies an altitude (m).  Z2 must not be
C  equal to Z1.

C  TV2 is a real variable that specifies the virtual temperature (K) at
C  Z2.

C  Z3 is a real variable that specifies an altitude (m).  Z3 may be
C  equal to Z1 or Z2.

C  RHO3 is a real variable that specifies the density (kg m**-3) at Z3.
C  Z1, TV1, Z2, TV2, Z3, and RHO3 define a density profile in a
C  constant-lapse-rate atmosphere.

C  Z_TOP_RAIN is a real variable that specifies the altitude (m) of the
C  top of the layer in which the hydrometeors are assumed to be liquid.

C  Z_BOT_SNOW is a real variable that specifies the altitude (m) of the
C  bottom of the layer in which the hydrometeors are assumed to be ice.
C  Z_BOT_SNOW must be greater than or equal to Z_TOP_RAIN.

C  Z is a real variable that specifies the altitude (m) at which an
C  echo-power-weighted hydrometeor terminal fall speed estimate is
C  desired.

C  DBZE is a real variable that specifies the measured radar
C  reflectivity factor (dBZ) for which an echo-power-weighted
C  hydrometeor terminal fall speed estimate is desired.

C  Output:

C  F is a real variable that returns an estimate of the
C  echo-power-weighted hydrometeor terminal fall speed.

C  SD_F is a real variable that returns the standard deviation of F.

      IMPLICIT NONE
      REAL DENCLR,TVCLR,DBX_TO_X
      REAL Z1,TV1,Z2,TV2,Z3,RHO3,Z_TOP_RAIN,Z_BOT_SNOW,Z,DBZE,F,SD_F,
     $Z_LIQ,Z_ICE,RHO,F_RAIN,F_SNOW,SD_F_RAIN,SD_F_SNOW,DEN_FACTOR,TV,
     $GAMMA

C  Calculate the air density.
      RHO=DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z)

C  Calculate the virtual temperature.
      TV=TVCLR(Z1,TV1,Z2,TV2,Z)

C  Convert the reflectivity factor from decibels for liquid and ice
C  scatterers.
      Z_LIQ=DBX_TO_X(DBZE)
      Z_ICE=DBX_TO_X(DBZE)*4.72

C  Calculate the modification for density.
      DEN_FACTOR=(1.+0.0023*(1.1-RHO/1.2040)*(293.15-TV))*
     $10.**(0.43*ALOG10(1.2040/RHO)-0.4*(ALOG10(1.2040/RHO))**2.5)

C  Calculate the terminal fall speed and its standard deviation for
C  solid hydrometeors.
      F_SNOW=DEN_FACTOR*0.817*Z_ICE**0.063
      IF(F_SNOW.LE.1.)THEN
         SD_F_SNOW=-0.5+1.*1.
      ELSE
         SD_F_SNOW=-0.5+1.*F_SNOW
      ENDIF

C  Calculate the terminal fall speed and its standard deviation for
C  liquid hydrometeors.
      GAMMA=(720.*0.08/(1.E-12*Z_LIQ))**(1./7.)
      F_RAIN=DEN_FACTOR*(9.65-10.30*(GAMMA/(GAMMA+6.))**7)
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
