      !    "<23-Jun-1994 19:28:02><UTC>"
      Subroutine Plt_Coverage (RangeKm,LuFile,
     > Olat,Olon,Itime_limits,LuMes)
c  Subroutine to plot the P-3 tail coverage:
c     only "Call Plot (x,y,ipen)" done here; other plot-related commands
c     should be done in calling routine (e.g., to set line width).
c
c  Read Lu LuFile for p3 tail coverage lines of form 
c    (F7.3,F9.3,4F6.1,1X,I4.4,2I2.2,"/",3I2.2)
c  and plot "tail beams" for start/stop YYYY,MM,DD,HH,MM,SS in
c  array Itime_Limits.

c  Input conditions:
c    LuFile: lu for file already opened and positioned for first
c        line of coverage (past header lines).
c    LuMes: lu for messages to user.
c    Itime_limits: array (dimensioned to 12) with start/stop times
c        YYYY,MM,DD,HH,MM,SS.
c    RangeKm: "range" in km for tail beam.  Often RangeKm will
c        be either about 38.4 (for 256 gates at 75m and next 128
c        gates at 150 m) or 76.8 km (for all gates) when the
c        variable range spacing was used.
c    Olat,Olon: origin latitude/longitude.
      
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      Implicit none

* parameters

      integer*4 LuFile, Lumes
     >     ,Itime_Limits(12)
      real*4  RangeKm ,Olat,Olon 


* locally defined parameters

      Integer*4 Itime(6)             ! Time read
     >  , KountRays                  ! count of rays read within time window
     > ,Iy,In,Id,Ih,Im,Is            !dates/time read
     > ,Iys, Ins, Ids, Ihs, Ims, Iss ! start time
     > , Iye, Ine, Ide, Ihe, Ime, Ise! end time
     > , Ipen                        ! =2 pen down, =3 pen up: for Plot calls
     > , Ierr                        ! for status returned from I/O calls
     > , i                           ! temp integer

      Logical*4 TimeOnBefore
      Real*4
     >  Plat,Plon,Azm        !lat,long,azmCorrected for roll
     >  ,GrAz,GrEl,Tilt      !ground-based azimuth("track")& elevation, tilt
     >  ,Ctr                 !factor to convert degrees to radians
     >  ,Dlat,Dlon           !difference in lat/lon
     >  ,DeltaX,DeltaY       !difference in X,Y
     >  ,TrigAng,Factor      !values used to help calculate X & Y
     >  ,X,Y                 !point for Plot call (values in kilometers)

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c scaling factor

      Ctr = 3.1415926535897/180.0


c  Define the start time for rays to read

      Iys = Itime_Limits(1) 
      Ins = Itime_Limits(2)
      Ids = Itime_Limits(3)
      Ihs = Itime_Limits(4)
      Ims = Itime_Limits(5)
      Iss = Itime_Limits(6)

c  Define the end time for rays to read

      Iye = Itime_Limits(7) 
      Ine = Itime_Limits(8)
      Ide = Itime_Limits(9)
      Ihe = Itime_Limits(10)
      Ime = Itime_Limits(11)
      Ise = Itime_Limits(12)



      Write (LuMes,100) Iys, Ins, Ids, Ihs, Ims, Iss, Iye, Ine, Ide,
     #                  Ihe, Ime, Ise
 100  Format (/'Plotting the tail coverage from ',i4,'/',i2.2,'/',
     #         i2.2,1x,3i2.2,' to ', i4,'/',i2.2,'/',i2.2, 1x,3i2.2/)

      KountRays = 0
      Ipen = 3

c  Read 1 second worth of data

 1    Read (LuFile,
     # '(F7.3,F9.3,4F6.1,1X,I4.4,2I2.2,"/",3I2.2)',
     #  Err=3,Iostat=Ierr,End=3)
     #  Plat,Plon,Azm,GrAz,GrEl,Tilt,
     #    (Itime(i),i=1,6)

c  Make sure the data is within the time window at the designated rate

      Iy = Itime(1) 
      In = Itime(2)
      Id = Itime(3)
      Ih = Itime(4)
      Im = Itime(5)
      Is = Itime(6)

      If (     TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iys,Ins,Ids,Ihs,Ims,Iss))
     #       Go To 1
      If (.not.TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iye,Ine,Ide,Ihe,Ime,Ise))
     #       Go To 2


  
  
      KountRays = KountRays + 1
      Dlat = Plat - Olat
      Dlon = Plon - Olon
      X = 111.19 * Cos(Plat*Ctr) * Dlon
      Y = 111.19 * Dlat
      TrigAng = (90.0 -GrAz) ! convert compass angle to trig angle

      Call Plot(x, y, Ipen)
      Ipen=3
      Factor = RangeKm * Cos(GrEl *Ctr)
      DeltaX = Cos (TrigAng*Ctr) * Factor
      DeltaY = Sin (TrigAng*Ctr) * Factor
      Call Plot(x+DeltaX, y+DeltaY, Ipen)
      Ipen=2
      Call Plot(x, y, Ipen)


      Go To 1

 2    Write (LuMes,'(" End time, Exiting p3-tail coverage plotting."
     # /" There were",I7," rays in time window."/)') KountRays

      Return

 3    if (ierr .ne. -1) Then
	  Write (LuMes,'("Error: ",i5," reading p3-tail"
     # ," coverage file.")') Ierr
      else
	  Write (LuMes,'("EOF in reading p3-tail coverage file.")') 
      end if
      Write (LuMes,'(" There were",I7," rays in time window."
     # /)') KountRays
      Return

      End !  subroutine Plt_Coverage ends
************************************************************************
