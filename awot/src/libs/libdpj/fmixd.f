      Function Fmixd (Dp,Press)

C   Returns MIXING RATIO IN                                (G/KG)
C   DP      DEWPOINT                                       (DEG C)
C   PRESS = PRESSURE                                       (MB)

      E = VAPOR(DP)
      FMIXD = 621.98*E/(PRESS-E)

      Return
      End
