      Function POTM(TEMP,PRESS,WMR)

C  Returns THE POTENTIAL TEMPERATURE OF MOIST AIR
C  TEMP  = TEMPERATURE                                    (DEG C)
C  PRESS = PRESSURE                                       (MB)
C  WMR   = MIXING RATIO                                   (G/KG)

      A    = 0.28544*(1.0-0.245*WMR*0.001)
      POTM = (TEMP+273.16)*(1000.0/PRESS)**A - 273.16

      Return
      End
