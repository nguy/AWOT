      Function Dewpt (T,Rh)

c  Function to return the dewpoint given the Temp (C) and
c    Relative Humidity (%)

      Q = Vapor(T) * Rh/100.0
      Enl = Alog(q)
      Dewpt = (243.5*Enl - 440.8)/(19.48-Enl)
      Return
      End
