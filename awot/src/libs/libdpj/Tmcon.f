      Function TMCON(TIME)  !Time: hhmmss->seconds

C Routine to convert time from hhmmss. format to seconds.

      IH = TIME/10000.0
      IM= (TIME-Float(IH)*10000.0)/100.0
      IS=TIME-Float(IH)*10000.0-Float(IM)*100.0
      TMCON=TIMZ(IH,IM,IS)

      Return
      End
