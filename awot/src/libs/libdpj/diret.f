      Function DIRET(U,V)

C   Returns DIRECTION (TOWARDS) GIVEN U & V COMPONENTS

       PI05 = 1.5707963
       PI15 = 4.7123889
       A    = PI05
       If (U.lt.0.0) A = PI15
       If (U.eq.0.0) U = .1E-32
       DIRET = 57.2958*(A-ATAN(V/U))
       Return
       End
