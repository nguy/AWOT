      Function POTT(TEMP,PRESS)

C   Returns THE POTENTIAL TEMPERATUE                       (DEG C)
C   TEMP  = TEMPERATURE                                    (DEG C)
C   PRESS = PRESSURE                                       (MB)

      POTT = (TEMP+273.16)*(1000.0/PRESS)**0.286714-273.16

      Return
      End
