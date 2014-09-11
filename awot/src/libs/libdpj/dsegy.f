      Function DSEGY(TEMP,RA,ALAT)
 
C**** THIS Function COMPUTES THE DRY STATIC ENERGY IN UNITS     (10**3 J/KG)
C**** DSE  = CP*T+PHI
C**** TEMP = AIR TEMPERATURE                                    (DEG.C)
C**** RA   = RADAR ALTITUDE                                     (M)
C**** ALAT = LATITUDE                                           (RADIAN)
 
      G     = 9.80616*(1.0 - 0.0026373*Cos(2.*ALAT))
      RE    =  G*648201.1446
      PHI   = (G*RE*RA)/(RE+RA)
      H     = 1004.0*(TEMP + 273.15)
      DSEGY =(H + PHI)*0.001
      Return
      End
