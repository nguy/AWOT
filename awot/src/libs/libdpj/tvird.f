      Function Tvird(Temp,Press,Dp)

C   Returns THE VIRTUAL TEMPERATURE                        (DEG C)
C   TEMP  = TEMPERATURE                                    (DEG C)
C   PRESS = PRESSURE                                       (MB)
C   DP    = DEWPOINT                                       (DEG C)

      W = Fmixd(Dp,Press)*0.001
      Tvird = (Temp+273.16)*((1.0+1.609*W)/(1.0+W)) - 273.16

      Return
      End
