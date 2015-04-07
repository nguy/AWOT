C  "<08-Oct-1993 21:54:35><UTC>"
**************************************************************************
       Logical*4 Function TimeOnBefore(
     >                  Iy1,Imon1,Id1,Ih1,Imin1,Isec1,
     >                  Iy2,Imon2,Id2,Ih2,Imin2,Isec2)! is time1<=time2?

c
c  No check is done for the inputs: year, month, day, hour, minute, and sec.
c  The year must be >= 1800.
c  TimeOnBefore returns .false.  if time1 is not on-or-before time2,
c               returns .true.  if time1 is on-or-before time2.
c
      Implicit None
      Integer*4 Iy1,Imon1,Id1,Ih1,Imin1,Isec1
      Integer*4 Iy2,Imon2,Id2,Ih2,Imin2,Isec2
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      TimeOnBefore = .true.   ! in case on-or-before
      if (Iy1 .lt. Iy2)then
	  return
      else if (Iy1 .gt. Iy2)then
          TimeOnBefore = .false.   ! not  on-or-before
	  return
      end if
      if (Imon1 .lt. Imon2)then
	  return
      else if (Imon1 .gt. Imon2)then
          TimeOnBefore = .false.   ! not  on-or-before
	  return
      end if
      if (Id1 .lt. Id2)then
	  return
      else if (Id1 .gt. Id2)then
          TimeOnBefore = .false.   ! not  on-or-before
	  return
      end if
      if (Ih1 .lt. Ih2)then
	  return
      else if (Ih1 .gt. Ih2)then
          TimeOnBefore = .false.   ! not  on-or-before
	  return
      end if
      if (Imin1 .lt. Imin2)then
	  return
      else if (Imin1 .gt. Imin2)then
          TimeOnBefore = .false.   ! not  on-or-before
	  return
      end if
      if (Isec1 .lt. Isec2)then
	  return
      else if (Isec1 .gt. Isec2)then
          TimeOnBefore = .false.   ! not  on-or-before
	  return
      end if
!The times are the same, so time1 <= time2
      TimeOnBefore = .false.   ! not  on-or-before
      return
      end ! Function TimeOnBefore ends
