Subroutine fchar(string)

  !  Routine to blank out non-charactes in a string

  Character String*(*)

  nch = Len(String)

  Do i = 1, nch
     If (String(i:i) .lt. ' ') String(i:i) = ' '
  End Do

End Subroutine fchar
