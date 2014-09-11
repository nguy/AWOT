      Function SfcP (T, Td, GA, P)
  
c  Routine to compute surface pressure by extrapolation of the d-value
c  to the surface using the standard atmospheric lapse rate

c   T  - Ambient Temperature (deg C)
c   TD - Dew Point Temperature (deg C)
c   P  - Static Pressure (mb)
c   GA - Geopotential altitude
c   SfcP - the estimated surface pressure

c  Calculate Temperature in K
      Temp = T + 273.16

c  What is the vapor pressure?
      W = Vp(Td)

c  What is the mixing ratio?
      Q = 0.622*W/(P-W)

c  What is the virtual temperature?
      Tv = Temp*((1.0+1.609*Q)/(1.0+Q))

c  What is the pressure altitude?
      Pa = PtoPA (P)

c  What would be the temperature at our height if we were in
c  at "standard atmosphere"?
      Tsa = 288.16 - 0.0065 * PA

c  Extrapolate the d-value (with a temperature correction) down to
c  z=0 to see what the surface pressure (SfcP) would be
      Pas =  PA-GA * TSA/TV 
      Sfcp = 1013.25*(1.0-PAS/44331.0)**5.25588

      Return
      End
