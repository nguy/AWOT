C  "<08-Oct-1993 21:54:33><UTC>"
************************************************************************
      Logical*4 Function OkTime(Iy,Imon,Id,Ih,Imin,Isec)!Check if time ok
c
c  Check is done for the input year, month, day, hour, minute, and sec.
c  The year must be >= 1800.
c  OkTime returns .false. if time is bad, .true. if it is ok.
c
      Implicit None
      Integer*4 Iy,Imon,Id,Ih,Imin,Isec
      Integer*4 Idays,Idaysmon
      Dimension Idays(12)
      Data Idays /31,28,31,30,31,30,31,31,30,31,30,31/
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      If ((Iy.lt.1800) .or. (Imon .lt.1).or.(Imon.gt.12))then
	  OkTime = .false.
	  Return 
      End If
  
      Idaysmon= Idays(Imon)
      If ((Imon.eq.2).and.(Mod(Iy,4).eq.0).and.
     > ((Mod(Iy,400).eq. 0).or.(Mod(Iy,100).ne.0)))Idaysmon=29! leap year
      if ((Id .lt.1) .or. (Id .gt. Idaysmon) .or.
     >    (Ih .lt.0) .or. (Ih .gt. 23) .or.
     >    (Imin .lt.0) .or. (Imin .gt. 59) .or.
     >    (Isec .lt.0) .or. (Isec .gt. 59))then
	  OkTime = .false.
	  Return 
      end if
      OkTime = .true.
      return
      end ! Logical*4 function OkTime ends
