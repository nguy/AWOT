      Function RHUMW(TEMP,PRESS,WMR)

C   Returns RELATIVE HUMIDTY                               (%)
C   TEMP  = TEMPERATURE                                    (DEG C)
C   PRESS = PRESSURE                                       (MB)
C   WMR   = MIXING RATIO                                   (G/KG)

      RHUMW = 100.0*WMR/STMIX(TEMP,PRESS)
      Return
      End
