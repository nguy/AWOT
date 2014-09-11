      Function Tvirr(Temp,Press,Rh)

C   TEMP  = TEMPERATURE                                    (DEG C)
C   PRESS = PRESSURE                                       (MB)
C   RH    = RELATIVE HUMIDITY                              (%)

      W = Fmixr(Temp,Press,Rh)*0.001
      Tvirr = (TEMP+273.16)*((1.0+1.609*W)/(1.0+W)) - 273.16

      Return
      End
