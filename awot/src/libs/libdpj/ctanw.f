      Function CTANW(TANW,RADW,ANGLE)

C Function to correct tangential wind for band crossing angle

      A=ANGLE*0.01745329
      CTANW=TANW*Cos(A) - RADW*SIN(A)
      Return
      End
