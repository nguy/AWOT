      Function FMIXR(TEMP,PRESS,RH)

C   Returns MIXING RATIO IN                                (G/KG)
C   TEMP  = TEMPERATURE                                    (DEG C)
C   PRESS = PRESSURE                                       (MB)
C   RH    = RELATIVE HUMIDITY                              (%)

      FMIXR = STMIX(TEMP,PRESS)*RH*0.01
      Return
      End
