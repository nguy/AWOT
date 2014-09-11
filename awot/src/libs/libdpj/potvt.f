      Function POTVT(TEMP,PRESS,WMR)

C   Returns THE POTENTIAL VIRTUAL TEMPERATURE              (DEG C)
C   TEMP  = TEMPERATURE                                    (DEG C)
C   PRESS = PRESSURE                                       (MB)
C   WMR   = MIXING RATIO                                   (G/KG)

      w = wmr*0.001
      TV = (TEMP+273.16)*((1.0+1.609*W)/(1.0+w)) - 273.16
      POTVT = POTT(TV,PRESS)

      Return
      End
