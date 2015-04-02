      Function ENERY(TEMP,DEWP,ADJT,PRES,PREF,RA,ALAT,ISEL)
 
C**** THIS Function COMPUTES THE DRY OR MOIST STATIC ENERGY IN UNITS
C****                                                           (1000*J/
C**** ENERGY = CP*T+PHI+L*W
C**** TEMP = AIR TEMPERATURE                                    (DEG.C)
C**** DEWP = DEWPOINT                                           (DEG.C)
C**** ADJT = TEMPERATURE ADJUSTMENT  TO THE REFERENCE LEVEL     (DEG.C)
C**** PRES = AIR PRESSURE                                       (MB)
C**** PREF = REFERENCE PRESSURE                                 (MB)
C**** RA   = RADAR ALTITUDE                                     (M)
C**** ALAT = LATITUDE                                           (RADIAN)
C**** ISEL = 0   DRY STATIC ENERGY ReturnED
C**** ISEL NON ZERO MOIST STATIC ENERGY ReturnED
C**** IF ADJT = 0.0 NO ADJUSTMENT IS MADE.
C**** WRITTEN BY KENNETH C. BELLE  NOAA/NHEML 02/25/80
C
      TP = TEMP+273.16+ADJT
      DP = DEWP+273.16+ADJT
      If ((Dp-Tp) .gt. 0.5 .and. ADJT .NE .0.0 .OR. (DP-TP) .gt. 0.5
     #  .and. ISEL .ne. 0) Go To 3
      AW = VP(DP-273.16)
      AW  = 0.622*AW /(PRES-AW )
      TV = TP*((1.0+1.609*AW)/(1.0+AW))
      G = 9.80616*(1.-0.0026373*Cos(2.*ALAT))
      ADRA = 0.0
      If (ADJT.eq.0.0) Go To 1

C**   GEOPOTENTIAL IS BEING ADJUSTED TO THE REFERENCE LEVEL

      ARA = 286.998*ALOG(PREF/PRES)
      ADRA = -(ARA*TV)/(G-ARA*0.5*0.005577)
    1 RH = RA + ADRA
      RE =  G*648201.1446
      PHI = (G*RE*RH)/(RE+RH)
      H = 1004.0*TP
      HL = 0.0
      If (ISEL.eq.0) Go To 2

C**   MOIST STATIC ENERGY IS REQUIRED SO THE LATENT HEAT IS BEING COMPUT

      HL = (2500.-2.274*(TP-273.16))*1000.*AW

    2 ENERY =(H+PHI+HL)*0.001
      Return

    3 ENERY = 1.E36
      Return
      End
