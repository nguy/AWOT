      Function WSEMR(TEMP,QRAT,RA,ALAT)
C
C**** THIS Function COMPUTES THE DRY STATIC ENERGY IN UNITS    (10**3 J/KG)
C**** WSE = CP*T+PHI+L*W
C**** TEMP = AIR TEMPERATURE                                    (DEG.C)
C**** QRAT = MIXING RATIO                                       (G/KG)
C**** PRES = AIR PRESSURE                                       (MB)
C**** RA   = RADAR ALTITUDE                                     (M)
C**** ALAT = LATITUDE                                           (RADIAN)
C
C
      If (QRAT.lt.0.) Go To 2
      AW    = QRAT*.001
      G     = 9.80616*(1.-0.0026373*Cos(2.*ALAT))
      RE    =  G*648201.1446
      PHI   = (G*RE*RA)/(RE+RA)
      H     = 1004.0*(TEMP + 273.15)
      HL = (2500.-2.274*TEMP)*1000.*AW
      WSEMR =(H+PHI+HL)*0.001
      Return
    2 WSEMR = 1.E36
      Return
      End
