      Function VAPOR(TT)

C   Returns VAPOUR PRESSURE   IF  TT = DEWPOINT            (MB)
C   Returns SATURATED VAPOUR PRESSURE  IF TT = TEMPERATURE (MB)
C   TT  = TEMPERATURE                                      (DEG C)

      T= TT+273.16

      If (t .lt. 0.0) Then
         vapor = 0.0
         Return
      End if
 
      VAPOR = 10.0**(22.5518-(2937.4/T)-4.9283*ALOG10(T))*10.
      Return
      End
