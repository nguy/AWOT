      Function RHUMD(TEMP,PRESS,DP)

C   Returns RELATIVE HUMIDITY                              (%)
C   TEMP  = TEMPERATURE                                    (DEG C)
C   PRESS = PRESSURE                                       (MB)
C   DP    = DEWPOINT                                       (DEG C)

      RHUMD = 100.*FMIXD(DP,PRESS)/STMIX(TEMP,PRESS)

      Return
      End
