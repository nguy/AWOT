
      Subroutine Change_Time(Iy,Imon,Id,Ih,Imin,Imin_Chan)!Domino time change
c
c  This is the same as Subroutine Change_TimeI2 except here all the
c     parameters passed are declared Integer*4.
c  This takes the change in current time in minutes (up to +/- 1440 minutes=
c     24 hours) and updates the year, month, day, hour, and minute.  This
c     works for the years 1800- (including when there are no leap years
c     in 1800, 1900, 2100, etc.).
c  No checking is done for the input year, month, day, hour, and minute.
c  If the minute change is greater than 1440, the routine will abort.
c
      Implicit None
      Integer*4 Iy,Imon,Id,Ih,Imin,Imin_Chan
      Integer*4 Idays,Idaysmon,Ihc,Idc
      Integer*4 LuWrite
      Parameter (LuWrite=6) ! Lu for abort message
      Dimension Idays(12)
      Data Idays /31,28,31,30,31,30,31,31,30,31,30,31/
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      If (abs(Imin_chan).gt.1440)Then
          Write (LuWrite,'("Abort in subroutine Change_time "
     > ,"  <08-Oct-1993 21:54:28><UTC>"
     1      /"because Abs(minute change)=",I7,"> 1440.")')Imin_chan
          Stop
      End If
  
      Imin = Imin + Imin_chan
      If((0.le.Imin).and.(Imin.le.59)) Return ! no more changes needed
      Ihc = Imin/60           ! hours to change
      Imin = Mod (Imin,60)
      If (Imin .lt. 0) Then
          Ihc = Ihc-1
          Imin = Imin + 60 ! put Imin in 0..59 range
      End If
      Ih = Ih + Ihc
      If ((0.le.Ih).and.(Ih.le.23)) Return ! no more changes needed
 
      Idc = Ih/24            ! days to change
      Ih = Mod(Ih,24)
      If (Ih.lt.0) Then
          Idc = Idc -1
          Ih = Ih +24 ! put Ih in 0..23 range
      End If
      Id = Id + Idc
      Idaysmon= Idays(Imon)
      If ((Imon.eq.2).and.(Mod(Iy,4).eq.0).and.
     > ((Mod(Iy,400).eq. 0).or.(Mod(Iy,100).ne.0)))Idaysmon=29! leap year
      If ((1.le.Id).and.(Id.le.Idaysmon)) Return ! no more changes needed
  
      If (Id .lt. 1) Then ! go to end of previous month
          Imon = Mod(Imon+10,12)+ 1! forward 11 = back 1
          Id = Idays(Imon)
          If ((Imon.eq.2).and.(Mod(Iy,4).eq.0).and.
     > ((Mod(Iy,400).eq. 0).or.(Mod(Iy,100).ne.0))) then !leap year
              Id=29
          Else If (Imon.eq.12) Then
              Iy=Iy-1
          End If
      Else
  
c  Here we need to go to the next month.
  
          Imon = Mod (Imon,12) + 1
          Id = 1
          If (Imon.eq.1) Iy = Iy + 1 ! Happy New Year!
      End If
      Return
      End ! subroutine Change_time ends
