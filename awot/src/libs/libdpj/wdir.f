      Function WDIR(U,V)

C  Returns THE WIND DIRECTION (FROM) GIVEN THE U & V COMPONENTS

       PI15 = 4.7123889
       PI05 = 1.5707963
       A    = PI15
       If (U.lt.0.0) A = PI05
       If (U.eq.0.0) U = .1E-32
       WDIR = 57.2958*(A-ATAN(V/U))
       Return
       End
