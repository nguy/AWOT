      !    "<23-Jun-1994 22:20:41><UTC>"
      Subroutine Plt_Map(LuFile,NameFile,
     >Clat, Clon, Xmin, Xmax, Ymin, Ymax,ScX,ScY, LuMes)
c  Subroutine to plot a map:
c     only "Call Plot (x,y,ipen)" done here; other plot-related commands
c     should be done in calling routine (e.g., to set line width).
c
c  Use Lu LuFile to open up and read the map file (NameFile) and
c  plot the map.  The file is closed before file exit.  
c  Messages will be written to Lu LuMes.
c  The map file is a byte stream (i.e., direct access) with multiple
c  entries (latitude,longitude,pen), where lat&lon are C float (assumed
c  here to be Real*4) and pen is C short-int (assumed here to be integer*2).
c  The byte-length per point should be BytesPerRec defined below.

c  Other input parameters
c  Clat, clon: bottom left corner latitude,longitude (origin)
c  Xmin, Xmax: X min/max (kilometers) of map to plot
c  Ymin, Ymax: Y min/max (kilometers) of map to plot
c    x=0,y=0 will correspond to clat,clon  
c  ScX, ScY: scale factors (normally = 1.0) to apply before plotting points.
c  
c Please Note: the following conditions will probably not plot correctly
c    (1) When two successive points are both outside the X-Y min/max,
c        then nothing will be plotted even though a portion of the
c        connecting line-segment may fall within the X-Y min/max.
c        This is usually not a serious problem if the data points are
c        close together.
c    (2) Any plot crossing 180 degrees west/east will probably only work
c        for the points with same longitude sign as Clon.
c        Suggested way to get around this is to make 2 calls to this 
c        routine, the 2nd one with Clon +/- 360.  E.g., if first call
c        is with Clon = 175, then 2nd call would be with 175-360=-185.
c        Or if first call is with Clon = -179, then 2nd call would be
c        with Clon = 181.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

      Implicit none

* parameters

      integer*4 LuFile, Lumes
      character NameFile*(*)
      Real*4 Clat,Clon,ScX,ScY,Xmin,Xmax,Ymin,Ymax

* local variables

      real*4
     >       Plat,Plon          ! most recent lat/long read
     >      ,Xkm1,Ykm1          ! point converted to Km coordinates
     >      ,Xkm2,Ykm2          ! most recent point converted to Km coordinates
     >      ,X,Y                ! point converted to Km coordinates scaled
      
      integer*4  KountPoints    ! number of points in map file read
     >         , KountPointsIn  ! number of points in map file read in window
     >         , Ipen           ! pen number
     >         , iout           ! help with clipping: =1 out of window,=0 in
     >         , ii             ! temp integer
     >         , ierr           ! temp integer for returning I/O error numbers
     >         , BytesPerRec    ! bytes for lat/long/pen in map file.
     >         , FirstNonChar   ! function to get first non-matching character
     >         , Ichk           !function to check where point is rel to min/max
     >         , FortranEOF     ! HP Fortran error code for reading past EOF
      Parameter (BytesPerRec=4+4+2)    ! bytes for lat/long/pen in map file.
      Parameter (FortranEOF=922) ! HP Fortran error code for reading past EOF
      integer*2  IpenRead       ! pen number read; integer*2 to save on space
      character file*100

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


      KountPoints = 0
      KountPointsIn = 0
      Ipen = 3

      ii = Max0(1,FirstNonChar(NameFile,' '))
      File = NameFile(ii:)  !use ii in case leading blanks
      Open (LuFile,Err=1000,File=File,Iostat=Ierr,Form='Unformatted',
     > Access='Direct',Recl=BytesPerRec,Status='Old')
c  Read 1 lat/long/pen triplet

      Read (LuFile,err=3,end=3,Iostat=ierr) Plat,Plon,IpenRead
      Ipen = IpenRead   ! convert to integer*4
      KountPoints = KountPoints + 1
  
      Iout = 1
      Call Convt_Km(Plat,Plon,Xkm1,Ykm1,Clat,Clon)
  
c  If first point is inside domain move pen with pen up
  
      If (Ichk(Xkm1,Ykm1,Xmin, Xmax, Ymin, Ymax) .eq. 0) Then
         KountPointsIn = KountPointsIn + 1
         Iout = 0
         x = Xkm1 * ScX
         y = Ykm1 * ScY
         Call Plot (x,y,3)
      End If
 
c  Read a coordinate from the disc file
  
    1 Read (LuFile,Err=3,End=3,Iostat=Ierr) Plat, Plon, IpenRead
      KountPoints = KountPoints + 1
      Ipen = IpenRead   ! convert to integer*4
  
c  Convert lat,lon to km
  
      Call Convt_Km(Plat,Plon,Xkm2,Ykm2,Clat,Clon)
  
c  Inside or outside of the domain?
  
      If (Ichk(Xkm2,Ykm2,Xmin, Xmax, Ymin, Ymax) .ne. 0) Go To 100
 
      KountPointsIn = KountPointsIn + 1
      If (Iout .eq. 1) Then
         If (Ipen .eq. 2) Call Clipit(Xkm1,Xkm2,Ykm1,Ykm2
     >    ,Xmin, Xmax, Ymin, Ymax,ScX,ScY)
         Xkm1 = Xkm2
         Ykm1 = Ykm2
         x = Xkm1 * ScX
         y = Ykm1 * ScY
         Call Plot (x,y,3)
         Iout = 0
         Go To 1
      End If
 
      x = Xkm2 * ScX
      y = Ykm2 * ScY
      Call Plot (x,y,Ipen)
      Xkm1 = Xkm2
      Ykm1 = Ykm2
      Iout = 0
      Go To 1
 
  100 If (Iout .eq. 1) Go To 101
      If (Ipen .eq. 2) Call Clipit(Xkm1,Xkm2,Ykm1,Ykm2
     >    ,Xmin, Xmax, Ymin, Ymax,ScX, ScY)
      Call Convt_Km(Plat,Plon,Xkm2,Ykm2,Clat,Clon)
      x = Xkm2 * ScX
      y = Ykm2 * ScY
      Call Plot (x,y,3)
      Iout = 1
  
  101 Xkm1 = Xkm2
      Ykm1 = Ykm2
      Go To 1

 3    if ((ierr .ne. -1)  .and.(ierr .ne. FortranEOF)) Then
	  Write (LuMes,'("Error: ",i5," reading map file",a)') Ierr
     >      ,NameFile
      else
	  Write (LuMes,'("EOF in reading map file",A)') NameFile
      end if
      Write (LuMes,'(" There were",I7," points read in routine"
     > ," Plt_Map,"/" with",I7," points in window.")') KountPoints
     > ,KountPointsIn
      Close (LuFile)
      Return
 1000 Write (LuMes,'(/"Error ",I5," opening map file ",a)')Ierr,
     >   NameFile
      Return

      End !  Subroutine Plt_Map ends
************************************************************************
