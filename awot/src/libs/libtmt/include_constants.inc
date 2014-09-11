C  include_constants.inc

C  Thomas Matejka NOAA/NSSL 18 November 1997

C  This INCLUDE file defines general constants.  Constants are in units
C  of kilograms-meters-seconds-Kelvin-moles-radians unless noted
C  otherwise.
C-----------------------------------------------------------------------
C  Mathematical constants.

C  Pi.
      REAL PI
      PARAMETER(PI=3.141592654)

C  Ratio of radians to degrees.
      REAL RADDEG
      PARAMETER(RADDEG=PI/180.)

C  Ratio of degrees to radians.
      REAL DEGRAD
      PARAMETER(DEGRAD=180./PI)

C  Gamma function at 1.5.
      REAL GAMMA_ONE_POINT_FIVE
      PARAMETER(GAMMA_ONE_POINT_FIVE=0.8862269)
C-----------------------------------------------------------------------
C  Physical and astronomical constants.

C  Speed of light in a vacuum.
      REAL C_LIGHT
      PARAMETER(C_LIGHT=2.99792E8)

C  Mean radius of the earth.
      REAL E_RAD
      PARAMETER(E_RAD=6371229.)

C  Angular velocity of the earth's rotation.
      REAL E_ROT
      PARAMETER(E_ROT=7.292E-5)

C  Gravitational acceleration at the earth's surface.
      REAL E_GRAV
      PARAMETER(E_GRAV=9.80616)
C-----------------------------------------------------------------------
C  Thermodynamic constants.

C  Molecular weight of dry air.
      REAL MW_DRY
      PARAMETER(MW_DRY=28.966E-3)

C  Molecular weight of water.
      REAL MW_H2O
      PARAMETER(MW_H2O=18.0160E-3)

C  Universal gas constant.
      REAL UNIV_GAS
      PARAMETER(UNIV_GAS=8.31436)

C  Specific gas constant of dry air.
      REAL R_DRY
      PARAMETER(R_DRY=UNIV_GAS/MW_DRY)

C  Specific gas constant of water vapor.
      REAL R_H2O
      PARAMETER(R_H2O=UNIV_GAS/MW_H2O)

C  Specific heat at constant pressure of dry air.
      REAL CP_DRY
      PARAMETER(CP_DRY=7.*UNIV_GAS/(2.*MW_DRY))

C  Specific heat at constant pressure of water vapor.
      REAL CP_VAP
      PARAMETER(CP_VAP=4.*UNIV_GAS/MW_H2O)

C  Specific heat at constant volume of dry air.
      REAL CV_DRY
      PARAMETER(CV_DRY=5.*UNIV_GAS/(2.*MW_DRY))

C  Specific heat at constant volume of water vapor.
      REAL CV_H2O
      PARAMETER(CV_H2O=3.*UNIV_GAS/MW_H2O)

C  Partial pressure of water at the triple point.
      REAL P_3PT
      PARAMETER(P_3PT=611.)

C  Temperature of water at the triple point.
      REAL T_3PT
      PARAMETER(T_3PT=273.16)

C  Melting/freeezing temperature of water in the presence of air.
      REAL T_FREEZE
      PARAMETER(T_FREEZE=273.15)

C  Density of ice.
      REAL RHO_ICE
      PARAMETER(RHO_ICE=917.)

C  Specific heat of ice.
      REAL C_ICE
      PARAMETER(C_ICE=2093.42)

C  Density of liquid water.
      REAL RHO_WAT
      PARAMETER(RHO_WAT=1000.)

C  Specific heat of liquid water.
      REAL C_WAT
      PARAMETER(C_WAT=4186.84)

C  Latent heat of sublimation of water at the triple point.
      REAL LH_SUB_3PT
      PARAMETER(LH_SUB_3PT=2.835E6)

C  Latent heat of evaporation of water at the triple point.
      REAL LH_VAP_3PT
      PARAMETER(LH_VAP_3PT=2.501E6)

C  Latent heat of melting of water.
      REAL LH_MELT
      PARAMETER(LH_MELT=LH_SUB_3PT-LH_VAP_3PT)

C  Reference pressure at which to define potential temperature.
      REAL P_REF
      PARAMETER(P_REF=100000.)
C-----------------------------------------------------------------------
C  Radar constants.

C  |K|**2 of ice.
      REAL K2_WAT
      PARAMETER(K2_WAT=0.93)

C  |K|**2 of liquid water.
      REAL K2_ICE
      PARAMETER(K2_ICE=0.197)

C  Rate of change of the natural logarithm of the index of refraction of
C  the atmosphere with respect to altitude.
      REAL DLNRIDZ
      PARAMETER(DLNRIDZ=-4.E-8)

C  Altitude of one atmosphere in the standard atmosphere.
      REAL,PARAMETER::Z0_SA=0.

C  Pressure of one atmosphere in the standard atmosphere.
      REAL,PARAMETER::P0_SA=101325.

C  Temperature at the altitude of one atmosphere in the standard
C  atmosphere.
      REAL,PARAMETER::T0_SA=288.15

C  Altitude of the tropopause in the standard atmosphere.
      REAL,PARAMETER::Z_TROP_SA=11000.

C  Temperature at the tropopause in the standard atmosphere.
      REAL,PARAMETER::T_TROP_SA=216.65

C  Lapse rate below the tropopause in the standard atmosphere.
      REAL,PARAMETER::GAMMA_TROP_SA=-0.0065

C  Lapse rate above the tropopause in the standard atmosphere.
      REAL,PARAMETER::GAMMA_STRAT_SA=0.

