      SUBROUTINE QP_FROM_DBZ(Z1,TV1,Z2,TV2,Z3,RHO3,Z_TOP_RAIN,
     $Z_BOT_SNOW,Z,DBZE,QP)

C  Thomas Matejka NOAA/NSSL 26 January 1997

C  This subroutine estimates the hydrometeor mixing ratio from the
C  reflectivity factor and altitude.

C  The relation for rain assumes an exponential distribution with the
C  Marshall-Palmer (1948) value for N0 = 0.08 cm**-4.

C  The relation for snow is based on results of Sekhon and Srivastava
C  (1970) (see Battan, 1973, p. 85).

C  The hydrometeor mixing ratio in the specified melting layer is that
C  of rain with the reflectivity factor reduced by a constant.

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

C  Z is a real variable that specifies the altitude (m) at which a
C  hydrometeor mixing ratio estimate is desired.

C  DBZE is a real variable that specifies the measured radar
C  reflectivity factor (dBZ) for which a hydrometeor mixing ratio
C  estimate is desired.

C  Output:

C  QP is a real variable that returns an estimate of the hydrometeor
C  mixing ratio.

      IMPLICIT NONE
      INCLUDE 'include_constants.inc'
      REAL,PARAMETER::DBZ_MELTING_OFFSET=8.
      REAL DENCLR,DBX_TO_X
      REAL Z1,TV1,Z2,TV2,Z3,RHO3,Z_TOP_RAIN,Z_BOT_SNOW,Z,DBZE,QP,Z_LIQ,
     $Z_ICE,RHO,GAMMA

C  Calculate the air density.
      RHO=DENCLR(Z1,TV1,Z2,TV2,Z3,RHO3,Z)

C  The hydrometeors are assumed to be ice.
      IF(Z.GE.Z_BOT_SNOW)THEN
         Z_ICE=DBX_TO_X(DBZE)*4.72
         QP=RHO_WAT*PI*0.000149*(1.E-12*Z_ICE)**0.3891/RHO

C  The hydrometeors are assumed to be liquid.
      ELSEIF(Z.LE.Z_TOP_RAIN)THEN
         Z_LIQ=DBX_TO_X(DBZE)
         GAMMA=(720.*0.08/(1.E-12*Z_LIQ))**(1./7.)
         QP=RHO_WAT*PI*0.08/(GAMMA**4*RHO)

C  The hydrometeors are assumed to be a mixture of ice and water.
C  Reduce the reflectivity factor and calculate as if liquid.
      ELSE
         Z_LIQ=DBX_TO_X(DBZE-DBZ_MELTING_OFFSET)
         GAMMA=(720.*0.08/(1.E-12*(Z_LIQ)))**(1./7.)
         QP=RHO_WAT*PI*0.08/(GAMMA**4*RHO)
      ENDIF

C  Done.
      RETURN
      END
