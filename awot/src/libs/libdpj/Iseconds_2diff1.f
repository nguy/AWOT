C  "<08-Oct-1993 21:54:30><UTC>"
*******************************************************************************
      Integer*4 Function Iseconds_2diff1(
     >                  Iy1,Imon1,Id1,Ih1,Imin1,Isec1,
     >                  Iy2,Imon2,Id2,Ih2,Imin2,Isec2)! time2-time1 in secs

c  Give the difference in seconds time2 -time1, where the difference must
c  be 0<=abs(diff)<=86400 (=24 hours), i.e., time1 must be within 24 hours of
c  time2.
c
c  No check is done for the inputs: year, month, day, hour, minute, and sec
c  to verify the times make sense.
c
      Implicit None
      Integer*4 Iy1,Imon1,Id1,Ih1,Imin1,Isec1
      Integer*4 Iy2,Imon2,Id2,Ih2,Imin2,Isec2
 
      Integer*4 IoneDay
      Parameter (IoneDay = 86400)! = 24 * 60 *60
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  This if test assumes there is at most 24 hours differences in the times.
      if (Iy2.gt.Iy1) then
	  Iseconds_2diff1 = IoneDay 
      else if (Iy2 .lt. Iy1) then
	  Iseconds_2diff1 = -IoneDay
      else if (Imon2 .gt. Imon1) then
	  Iseconds_2diff1 = IoneDay
      else if (Imon2 .lt. Imon1) then
	  Iseconds_2diff1 = -IoneDay
      else if (Id2 .gt. Id1) then
	  Iseconds_2diff1 = IoneDay
      else if (Id2 .lt. Id1) then
	  Iseconds_2diff1 = -IoneDay
      else
	  Iseconds_2diff1 = 0
      end if
      Iseconds_2diff1 = Iseconds_2diff1 + 3600*(Ih2-Ih1)+
     >   60*(Imin2-Imin1) + Isec2 -Isec1
      return
      end ! Integer*4 Function Iseconds_2diff1 ends
