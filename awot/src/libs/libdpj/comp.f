      Subroutine COMP(U,V,WD,WS)

c Routine to calculate u,v components from the wind direction [deg]
c and speed [m/s]

C WD WITH RESPECT TO NORTH

      U = Ws * Cos((270.0-Wd)*.0174532925)
      V = Ws * Sin((270.0-Wd)*.0174532925)

      Return
      End
