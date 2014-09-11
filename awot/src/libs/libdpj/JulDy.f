C  "<08-Oct-1993 21:54:31><UTC>"
**************************************************************
      Integer*4 Function Juldy(Imon,Id,Iy)! get Julian day

C*** COMPUTES JULIAN DAY FROM THE DATE (IMON,IDAY,IYEAR)
C*** ----------------------------------(  09, 27 , 1980)
c  Little checking is done for the input year, month, day,
c   except the year must be >= 1800.
c
      Implicit None
      Integer*4 Iy,Imon,Id
      Integer*4 ITD
      Dimension ITD(12)
      DATA ITD/31,59,90,120,151,181,212,243,273,304,334,365/
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      If ((Iy.lt.1800) .or. (Imon .lt.1).or.(Imon.gt.12))then
	  JulDy = -999
	  Return 
      End If
      JulDy = Id
      if (Imon .eq. 1) return ! for January that's all
      JulDy = JulDy + ITD(Imon-1)  ! add in previous months' days
      If ((Imon.ge.3).and.(Mod(Iy,4).eq.0).and.
     > ((Mod(Iy,400).eq. 0).or.(Mod(Iy,100).ne.0)))
     >  JulDy = JulDy + 1  ! add in for March or later in leap year
      Return
      End ! Integer*4 Function JulDy ends
