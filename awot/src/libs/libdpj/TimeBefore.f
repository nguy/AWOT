C  "<08-Oct-1993 21:54:35><UTC>"
**************************************************************************
       Logical*4 Function TimeBefore(
     >                  Iy1,Imon1,Id1,Ih1,Imin1,Isec1,
     >                  Iy2,Imon2,Id2,Ih2,Imin2,Isec2)! is time1< time2?

c
c  No check is done for the inputs: year, month, day, hour, minute, and sec.
c  The year must be >= 1800.
c  TimeBefore returns .false.  if time1 is after or equal to time2,
c               returns .true.   if time1 is before time2.
c
      Implicit None
      Integer*4 Iy1,Imon1,Id1,Ih1,Imin1,Isec1
      Integer*4 Iy2,Imon2,Id2,Ih2,Imin2,Isec2
c---------------------------------------------------------------

      TimeBefore = .true.        ! in case before

      If (Iy1 .lt. Iy2) Then
	  Return
      Else If (Iy1 .gt. Iy2) Then
          TimeBefore = .false.   ! not before
	  Return
      End If

      If (Imon1 .lt. Imon2) Then
	  Return
      Else If (Imon1 .gt. Imon2) Then
          TimeBefore = .false.   ! not before
	  Return
      End If

      If (Id1 .lt. Id2) Then
	  Return
      Else If (Id1 .gt. Id2) Then
          TimeBefore = .false.   ! not before
	  Return
      End If

      If (Ih1 .lt. Ih2) Then
	  Return
      Else If (Ih1 .gt. Ih2) Then
          TimeBefore = .false.   ! not before
	  Return
      End If

      If (Imin1 .lt. Imin2) Then
	  Return
      Else If (Imin1 .gt. Imin2) Then
          TimeBefore = .false.   ! not before
	  Return
      End If

      If (Isec1 .lt. Isec2) Then
	  Return
      Else If (Isec1 .gt. Isec2) Then
          TimeBefore = .false.   ! not before
	  Return
      End If

!The times are the same, so time1 = time2

      TimeBefore = .False.   ! not  before

      Return
      End
