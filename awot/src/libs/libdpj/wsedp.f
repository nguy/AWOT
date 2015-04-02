      Function WSEDP(TEMP,DEWP,ADJT,PRES,PREF,RA,ALAT)
C
C**** THIS Function COMPUTES THE DRY STATIC ENERGY IN UNITS    (10**3 J/KG)
C**** WSE  = CP*T+PHI+L*W
C**** TEMP = AIR TEMPERATURE                                    (DEG.C)
C**** DEWP = DEWPOINT                                           (DEG.C)
C**** ADJT = TEMPERATURE ADJUSTMENT  TO THE REFERENCE LEVEL     (DEG.C)
C**** PRES = AIR PRESSURE                                       (MB)
C**** PREF = REFERENCE PRESSURE                                 (MB)
C**** RA   = RADAR ALTITUDE                                     (M)
C**** ALAT = LATITUDE                                           (RADIAN)
C**** IF ADJT = 0.0 NO ADJUSTMENT IS MADE.
      TP    = TEMP + ADJT
      DP    = DEWP + ADJT
      If ((DP-TP).gt.0.5.and.ADJT.ne.0.0) Go To 2
      AW    = VAPOR(DP)
      AW    = 0.622*AW /(PRES-AW )
      TV    = (TP+273.15)*((1.0+1.609*AW)/(1.0+AW))
      G     = 9.80616*(1.-0.0026373*Cos(2.*ALAT))
      ADRA  = 0.0
      If (ADJT.eq.0.0) Go To 1
C
C**   GEOPOTENTIAL IS BEING ADJUSTED TO THE REFERENCE LEVEL
C
      ARA   = 286.998*ALOG(PREF/PRES)
      ADRA  = -(ARA*TV)/(G-ARA*0.5*0.005577)
    1 RH    = RA + ADRA
      RE    =  G*648201.1446
      PHI   = (G*RE*RH)/(RE+RH)
      H     = 1004.0*(TP + 273.15)
      HL    = (2500.-2.274*TP)*1000.*AW
      WSEDP =(H+PHI+HL)*0.001
      Return
    2 WSEDP = 1.E36
      Return
      End
