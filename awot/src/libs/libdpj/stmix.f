      Function STMIX(TEMP,PRESS)

C   Returns SATURATED MIXING RATIO                         (G/KG)
C   TEMP  = TEMPERATURE                                    (DEG C)
C   PRESS = PRESSURE                                       (MB)

      ES  = VAPOR(TEMP)
      STMIX = 621.98*ES/(PRESS-ES)
      Return
      End
