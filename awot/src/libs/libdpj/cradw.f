      Function CRADW(TANW,RADW,ANGLE)

C Function to correct radial wind for band crossing angle

      A=ANGLE*0.01745329
      CRADW=TANW*SIN(A) + RADW*Cos(A)
      Return
      End
