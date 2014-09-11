      Function SPHUM(DP,PRESS)

C   Returns SPECIFIC HUMIDTY IN                            (G/KG)
C   DP      DEWPOINT                                       (DEG C)
C   PRESS = PRESSURE                                       (MB)

      E = VAPOR(DP)
      SPHUM = 621.98*E/(PRESS-0.37803*E)
      Return
      End
