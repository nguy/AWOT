      Function Atn2(y,x)

c  Routine to compute an arc-tangent from two arguements

      If (x .lt. 0.0001 .and. x .gt. -0.0001) Go To 2  ! x=0 ?
      If (y .ge. 0.0) Go To 7                          ! y>0?

      Atn2 = 6.283185308 + Atan2(y,x)
      Return

 7    Atn2 = Atan2(y,x)
      Return

 2    If (y) 4,3,5

 3    Atn2 = 3.141592654
      Return

 4    Atn2 = 4.712388981
      Return

 5    Atn2 = 1.570796327
      Return
      End
