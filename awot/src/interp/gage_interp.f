      Program gage_interp

c  Routine to interpolate rain gage data to a Cartesian grid

      Common /Block/ Imax, Jmax, Sx, Sy, Xmax, Ymax, Clat, Clon,
     #               Iys, Ins, Ids, Ihs, Ims, Iss,
     #               Iye, Ine, Ide, Ihe, Ime, Ise

      PARAMETER (Islim = 100)  ! Max number of gages for the interpolation
      REAL X(ISLIM), Y(ISLIM), Z(ISLIM)
      Integer FirstChar

      Character*100 File, FileOut, Gages, Line, RainDir, ESRL_Gages*100
      Character StName*36, Sta*3

      Logical*4 TimeOnBefore
      Data Ctr /0.01745329/

C    Open the data file containing the parameter information.

      RainDir = '/Users/davej/programs/hmt/Rainfall/'
      ESRL_Gages ='    '

      n = IargC()
      
      If (n .eq. 0) Then
         Write (6,'("Enter prm file name:",$)')
         Read (5,'(a)') File
      Else
         Call GetArg(1, File)
      End If
 
      nch = LenTrim(File)
      File = File(1:nch)

C  Open the parameter file
  
      Open (10,Err=1020,File=File,Iostat=Ierr,Status='Old')

      Read (10,'(a)') Line
      Call CaseFold(Line)
      nch = LenTrim(Line)

      If (Line(1:11) .ne. 'GAGE_INTERP') Then
         Write (6,'("wrong parameter file - line 1 not gage_interp ",a)')
     $        Line(1:nch)
         Stop 1
      End If

c  Read the grid corner point and dimensions of the output grid
c     Clat,Clon - grid lower left corner coordinate
c     Imax,Jmax - grid size in points
c     Sx, Sy    - distance (km) between grid points

      Read (10,*) Clat, Clon, Imax, Jmax, Sx, Sy

c  Starting and ending times

      Read (10,*) Isymd, Ishms, Ieymd, Iehms

c  Grid size in km

      Xmax = Float(Imax) * Sx
      Ymax = Float(Jmax) * Sy

      Write (6,'("Grid Size (pts,km): ",2i4,1x,2F9.3)') Imax, Jmax, Xmax
     $     , Ymax 

      Call Tcnvt(Isymd, Iys, Ins, Ids)
      Call Tcnvt(Ishms, Ihs, Ims, Iss)
      Call Tcnvt(Ieymd, Iye, Ine, Ide)
      Call Tcnvt(Iehms, Ihe, Ime, Ise)

      Write (6,'("Starting Time: ",3i2.2,1x,3i2.2,1x,
     #           "Ending Time:   ",3i2.2,1x,3i2.2," UTC")')
     #     Iys, Ins, Ids, Ihs, Ims, Iss,
     #     Iye, Ine, Ide, Ihe, Ime, Ise

      Read (10,'(a)') FileOut
      Read (10,'(a)') Gages
      Read (10,'(a)',End=340) ESRL_Gages
 340  Close(10)

      nch = FirstChar(FileOut,' ') - 1
      FileOut = FileOut(1:nch)
      Write (6,'(/"Output file: ",a)') FileOut

      nch = FirstChar(Gages,' ') - 1
      Gages = Gages(1:nch)
      Write (6,'(/"Rain gages file: ",a)') Gages

      File = Gages
      Open (15, Err=1020, File=File, Iostat=Ierr, Status='Old')
      Read (15,*)

      Nsta = 0
      nchR = LenTrim(RainDir)

 350  Read (15,'(a)', End=88) Line

      Read (Line,'(a,a,10x,i7,2x,F8.3,1x,f9.3,4x,F6.1,2x,F6.2)')
     #     Stname, Sta, Ialt, Xlat, Xlon

      If (Stname(1:5) .eq. 'xxxxx') Go To 88

c  Open the rainfall file

      File = RainDir(1:nchR) // Sta // '.dat'
      
      Open (11, Err=1020, File=File, Iostat=Ierr, Status='Old')

c  Skip the first 3 lines

      Read (11,*)
      Read (11,*)
      Read (11,*)

c  Read the accumulation between the two time periods

 500  Read (11,1001,End=99, Iostat=Ierr) Iy, In, Id, Ih, Im
     $     , Rain
 1001 Format (I4,I2,I2,1x,I2,I2,1x,F6.0)

c     The raingage data from HMT is usually in local time (i.e., PST),
C     need to convert the times from the files to UTC to match the radar
C     data (and all the other data sets collected in UTC)
C     The conversion from PST to UTC is to ADD 8 hours to PST

      Is = 0
      Call Chtime(Iy, In, Id, Ih, Im, Is)

c  Find the first record past the start time

      If (TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iys+2000,Ins,Ids,Ihs,Ims,Iss)
     $     )Go To 500

      Rain1 = Rain

 501  Read (11,1001,End=99, Iostat=Ierr) Iy, In, Id, Ih, Im
     $     , Rain

c  Find the first record after the stop time

      Call Chtime(Iy, In, Id, Ih, Im, Is)

      If (TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iye+2000,Ine,Ide,Ihe,Ime,Ise)) 
     $     Go To 501

c  How much rain is the difference between the two amounts

      RainT = Rain - Rain1
      If (RainT .lt. 0.0) RainT = 0.0

      RainT = RainT * 25.4    ! Convert to mm from inches

      Ypos = (Xlat - Clat) * 111.19
      Xpos = (Xlon - Clon) * 111.19 * Cos(Xlat*Ctr)
      Nsta = Nsta + 1

      If (Nsta .gt. Islim) Then
         Write (6,'("Too many raingage stations",2i5)') Nsta, Islim
         Stop 999
      End If

c  x,y,and z are the input (irregularly spaced) gage data

      x(Nsta) = Xpos
      y(Nsta) = Ypos
      z(Nsta) = RainT

      Write (6,'(i2,1x,a,1x,"R1 [in]= ",F6.2," R2 [in]= ",F6.2
     $     ," Tot [mm]=",F7.2," x,y: ",2F9.3)') Nsta, Sta, Rain1, Rain,
     $     RainT, Xpos,Ypos 

      Close (11)
      Go To 350

c  program flow ends up here if now start or end time was found

 99   Write (6,'(/"No start or end time found for station ",a)') Sta
      Go To 350

 88   Close(15)

c  Are there any ESRL gages?

      nchE = FirstChar(ESRL_Gages,' ') - 1
      Write (6,'(/"ESRL Gage File: ",a)') ESRL_Gages(1:nchE)

      If (nchE .gt. 0) Then
         File = ESRL_Gages(1:nchE)
         Open (14, File=File, Iostat=Ierr, Status='Old')

 80      Read (14,'(a)') File
         If (File(1:4) .eq. 'xxxx') Go To 87

         Write (6,'("ESRL gage: ",a)') File
         Open (13, File=File, Err=1020, Iostat=Ierr, Status='Old')
         Read (13,'(a)') Line
         Read (13,'(a)') Line
         Write (Sta,'(a3)') Line(7:9)

         Read (13,'(a)') Line
         Read (Line,'(6x,F12.0)') Xlat
         Read (13,'(a)') Line
         Read (Line,'(6x,F12.0)') Xlon
         Read (13,'(a)') Line
         Read (13,'(a)') Line

         Ypos = (Xlat - Clat) * 111.19
         Xpos = (Xlon - Clon) * 111.19 * Cos(Xlat*Ctr)
         Nsta = Nsta + 1
         RainT = 0.0

 70      Read (13,'(5x,i3,4x,i2,i2,3x,F12.0)',End=71) JD, Ih, Im, Rain

         Iyear = 2005
         If (JD .lt. 350) Iyear = 2006

         Jday_Start = JulDate(Iyear, 1, 1)
         Jdate = Jday_Start + JD - 1
         Call Jdate2Day(JDate, Iy, In, Id)

c  Find the first record past the start time

         If (TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iys+2000,Ins,Ids,Ihs,Ims,Iss
     $        )) Go To 70

         If (.not.TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iye+2000,Ine,Ide,Ihe
     $        ,Ime,Ise)) Go To 71

         RainT = Rain + RainT
         Go To 70

 71      x(Nsta) = Xpos
         y(Nsta) = Ypos
         z(Nsta) = RainT
         
         Write (6,'(i3,1x,a,1x,
     $        " Tot [mm]=",F7.2," x,y: ",2F9.3)') Nsta, Sta,
     $        RainT, Xpos,Ypos

         Go To 80          
      End If

 87   Call InterpIt(Nsta, X, Y, Z, FileOut) 

      Call Exit

 1020 Write (6,'("Error on file:",a, "Err=",i5)') File, Ierr
      Stop 1020
      End

      Subroutine InterpIt(Nsta, X, Y, Z, FileOut)

      Common /Block/ Imax, Jmax, Sx, Sy, Xmax, Ymax, Clat, Clon,
     #               Iys, Ins, Ids, Ihs, Ims, Iss,
     #               Iye, Ine, Ide, Ihe, Ime, Ise

c     Routine to interpolate irregularly spaced gages to a Cartesian
C     grid using the NCAR "natgrid" routine

      Dimension XI(Imax), YI(Jmax), ZI(Imax,Jmax)
      Dimension X(*), Y(*), Z(*)

      Character Time_Current*24, Ident*8, Project*16
      Character FileOut*100

C  Define the output grid for interpolation

      Xmin =  0.0
      Ymin = 0.0
      XINC = (XMAX-XMIN)/Float(Imax-1)

      DO I = 1, Imax
         XI(I) = XMIN + Float(I-1) * XINC
      End Do

      YINC = (YMAX-YMIN)/Float(Jmax-1)

      DO J = 1, Jmax
         YI(J) = YMIN + Float(J-1) * YINC
      End Do

C  Set the flag for using estimated gradients.

      CALL NNSETI ('IGR',1)

c  Set the flag to not allow negative numbers

      CALL NNSETI('NON - control of negative values',1)

C  Do the gridding.

      CALL NATGRIDS (Nsta, X, Y, Z, Imax, Jmax, XI, YI, ZI, IER)

      IF (IER .NE. 0) THEN
         WRITE (6,510) IER
 510     FORMAT('Error return from NATGRIDS = ',I3)
         Stop 510
      END IF

c Write an output file

      Open (15, Err=1020, File=FileOut, Iostat=Ierr,
     #     Form='UNFORMATTED')

      Call FDATE(Time_Current)
      Ident = 'gage_amt'
      Project = 'HMT-05/06       '

      Write (15) Ident, Time_Current, Clat, Clon, 
     #     Sx, Sy, Imax, Jmax, Project, 
     #     Iys, Ins, Ids, Ihs, Ims, Iss, 
     #     Iye, Ine, Ide, Ihe, Ime, Ise
      
      Write (6,'(/"Writing output file: ",a)') FileOut
      
      Do j = 1, Jmax
         Write (15) (Zi(i,j),i=1,Imax)
      End Do
      
      Close (15)

      Return

 1020 Write (6,'("Error on file:",a, "Err=",i5)') FileOut, Ierr
      Return

      End

      Subroutine Chtime (Iy, In, Id, Ih, Im, Is)

c  Routine to add 8 hours to the time

      Jd = Juldate(Iy, In, Id)
      Tm = Timz(Ih, Im, Is)
      Tm = Tm + 28800.0           ! 8 hours in seconds
      Iday_Offset = 0

      If (Tm .ge. 86400.0) Then
         Tm = Tm - 86400.0        ! day has flipped
         Iday_Offset = 1
      End If

      Call Ctme(Tm, Ih, Im, Is)
      Jd = Jd + Iday_Offset

      Call Jdate2Day(Jd, Iy, In, Id)
      Return
      End
