      Subroutine Ctme (T,Ih,Im,Is)

c  Routine to convert time from seconds to hour,minutes,seconds

      Ih = T/3600.0                  ! Hours
      Im = (T-Float(Ih)*3600.0)/60.0 ! Minutes
      Is = T-Float(Ih)*3600.0-Float(Im)*60.0
      Return

      End
