      Program grnd_interp

!  Program to interpolate 88D radar date to 3-D grids
!  David Jorgensen NOAA/NSSL April 2006

!  All printed output is directed to a .log file
!  Radar data is input as Universal Format (UF)

        ! The UF file is put on the command line, e.g., grnd_interp prm_file UF_file
!  The parameters used in the program:

!  Namof- name of the disc file that will contain output Cartesian data
!  Namdf- name of the disc file that contains the UF 88D data
!  Imax - maximum points in X direction (East-West)
!  Jmax - maximum points in Y direction (North-South)
!  Kmax - maximum points is Z direction (vertical)
!  Sx   - distance between points (km) in X direction
!  Sy   - distance between points (km) in Y direction
!  Sz   - distance between points (km) in Z direction
!  Z0   - height above surface of 1st plane (km)
!         (Note: Definition of "surface" is defined as that
!         of mean sea level

!  Nmosm- flag for type of grid reference desired
!       - =0 with respect to the ground; in this case,
!          OLAT - latitude of lower left hand corner
!          OLON - longitude of lower left hand corner
!       - =1 with respect to the moving grid center; in that case,
!          Olat = Olat + Sv*(Time-T_init)/111190.0
!          OLon = Olon + Su*(Time-T_init)/(111190.0*Cos(Olat))
!  IHS, IMS, ISS - starting time of 1st leg
!  IHE, IME, ISE - ending time of 1st leg

!  Summary of disc lu's used:
!       10 - parameter file
!       60 - log file
!       50 - input UF file
!       95 - OUTPUT (Cartesian) file

      Integer :: FirstChar
      Integer, Parameter :: MaxG = 100

      Real, Dimension(MaxG) :: Vel, Ref
      Logical :: OKTime, EOF

      Common /One/ Imax, Jmax, Kmax
      Common /Two/ Sx, Sy, Sz, Z0, Olat, Olon, Nmosm, Su, Sv, Xniq, Alt, DZ_Name, Vr_Name
      Common /Five/ Itime_Limits(12), Init_Time(6), Isday
      Common /Write/ LM
      Common /Interp/ Hinf, Vint, Vslp

      Character(Len=100) :: Line, Ifile, PrmFile, Namdf, Namof, LogFile
      Character Time_Current*24, UserN*8, Rname*8, Project*16
      Character(Len=2) DZ_Name, Vr_Name
      Character(Len=3) Itype, Type

      nparms = IargC()
      
      If (nparms .lt. 1) Then
         Write (6,'("Enter prm file name:",$)')
         Read (5,'(a)') PrmFile
      Else
         Call GetArg(1,PrmFile)
      End If

      nchp = LenTrim(PrmFile)
      Ifile = PrmFile(1:nchp)

!    Open the parameter file

      Open (10, File=Ifile, Err=1010, Iostat=Ierr, Status='Old')
      
!  Read parameters from the disc file
!     The first line is the indicator for the type of .prm file (grnd_interp.prm)

      Read (10,'(a)') Line
      Call CaseFold(Line)
      nch = FirstChar(Line,' ') - 1

      If (Line(1:nch) .ne. 'GRND_INTERP') Then
         Write (6,'(/"  Improper type of .prm file - must start with ","Grnd_interp ",a/)') Line(1:nch)
         Stop 'Grnd_Interp'
      End If

!  Project name 

      Read (10,'(a)') Project

!  Get the names of the UF input file (Namdf) & output gridded file (Namof)

      Read (10,'(a)') Namdf
      Read (10,'(a)') Namof

!     Open the UF data file (Namdf)

      nchd = FirstChar(Namdf,' ') - 1
      Ifile = Namdf(1:nchd)

      Open (50,File=Ifile, Form='UNFORMATTED', Err=1010, Iostat=Ierr, Status='Old', Convert='Big_Endian')

!     Name the log file the same as the .prm file but replace the suffix with .log

      If (LM .ne. 6) Then
         nchl = FirstChar(PrmFile,' ') - 4
         LogFile = PrmFile(1:nchl) // 'log'
         Ifile = LogFile
         Open (LM, File=Ifile, Iostat=Ierr, Err=1010)
      End If

!     Open the output file

      ncho = FirstChar(Namof,' ') - 1
      Ifile = Namof(1:ncho)

      Open (95, File=Ifile, Iostat=Ierr, Err=1010, Form='Unformatted')  
      Write (LM,'(/"Output file name: ",a)') Namof(1:ncho)

!--------------------------------------------------------------
!  Write out info to the .log file, including the user name
!--------------------------------------------------------------

      UserN = 'xxxxxxxx'
      Call GetLog(UserN)

      Write (LM,'("Hybrid-Scan Radar Interpolation program - .prm file:",a, &
                  " for user: ",a)') Prmfile(1:nchp), UserN

!  The grid parameters

      Read (10,*) Imax, Jmax, Kmax, Sx, Sy, Sz, Z0
      Read (10,*) Nmosm, Olat, Olon, Iynd, Ihms, Su, Sv

      Call Tcnvt(Iynd, Iyi, Ini, Idi)            ! Convert from yymmdd to year, month, day
      Call Tcnvt(Ihms, Ihi, Imi, Isi)            ! Convert from yymmdd to year, month, day

!     Itype is the 3 character ID of the radar (e.g., 'OUN' for KOUN)
!     The altitude (MSL) of the radar. Need this because sometimes the
!     altitude from the UF file is weird, especially for the SMART-R.
!     DZ_Name is the name of the reflectivity field in the UF file (A2)
!     Vr_Name is the name of the radial velocity field in the UF file (A2)
!     Hint, Vint and Vslp are the cressman weighting intervals (horizontal, vertical, and slope)

      Read (10,*) Itype, Alt, DZ_Name, Vr_Name, Hinf, Vint, Vslp

!     Start time, End time

      Read (10,*) Iynds, Ihmss, Iynde, Ihmse

      Call Tcnvt(Iynds, Iys, Ins, Ids)     ! Convert from yymmdd to year, month, day
      Call Tcnvt(Ihmss, Ihs, Ims, Iss)     ! Convert from yymmdd to year, month, day
      Call Tcnvt(Iynde, Iye, Ine, Ide)     ! Convert from yymmdd to year, month, day
      Call Tcnvt(Ihmse, Ihe, Ime, Ise)     ! Convert from yymmdd to year, month, day

      Close(10)

      Write (LM,'(/"Input Variables:"/&
          " Project Name:",1x,a/      &
          " Imax:",1x,i4/             &
          " Jmax:",1x,i4/             &
          " Kmax:",1x,i4/             &
          " Sx:  ",1x,F6.3/           &
          " Sy:  ",1x,F6.3/           &
          " Sz:  ",1x,F6.3/           &
          " Z0:  ",1x,F6.3/           &
          " Nmosm:",1x,I2/            &
          " Olat: ",1x,F7.3/          &
          " Olon: ",1x,F8.3/          &
          " Init Time: ",1x,i6.6,1x,i6.6/ &
          " Su:   ",1x,F6.2/          &
          " Sv:   ",1x,F6.2/          &
          " Alt:  ",1x,F6.2/          &
          " DZ Field in UF file=",1x,a/             &
          " Vr field in UF file=",1x,a/             &
          " Hinf, Vint, Vslp",1x,3F7.2)')&
          Project, Imax, Jmax, Kmax, Sx, Sy, Sz, Z0, Nmosm, Olat, Olon, Iynd, Ihms, Su, Sv, Alt, DZ_Name, Vr_Name, Hinf, Vint, Vslp

!  Check the grid motion init time

      If (Nmosm .eq. 1) Then
         If (.not. OKTime(Iyi+2000,Ini,Idi,Ihi,Imi,Isi)) Then
            Write (LM,'("Time: ",i2.2,"/"i2.2,"/",i2.2,1x,3i2.2," is not OK")') Iyi, Ini, Idi, Ihi, Imi, Isi
            Stop
         End If
      End If

      Init_Time(1) = Iyi + 2000
      Init_Time(2) = Ini
      Init_Time(3) = Idi
      Init_Time(4) = Ihi
      Init_Time(5) = Imi
      Init_Time(6) = Isi

!  Are the start times OK?

      If (.not. OKTime(Iys+2000,Ins,Ids,Ihs,Ims,Iss)) Then
         Write (LM,'("Time: ",i2.2,"/"i2.2,"/",i2.2,1x,3i2.2," is not OK")') Iys,Ins,Ids,Ihs,Ims,Iss
         Stop
      End If

!  Are the end times OK?

      If (.not. OKTime(Iye+2000,Ine,Ide,Ihe,Ime,Ise)) Then
         Write (LM,'("Time: ",i2.2,"/"i2.2,"/",i2.2,1x,3i2.2," is not OK")') Iye,Ine,Ide,Ihe,Ime,Ise
         Stop
      End If

      Itime_Limits(1) = Iys + 2000
      Itime_Limits(2) = Ins
      Itime_Limits(3) = Ids
      Itime_Limits(4) = Ihs
      Itime_Limits(5) = Ims
      Itime_Limits(6) = Iss
      Itime_Limits(7) = Iye + 2000
      Itime_Limits(8) = Ine
      Itime_Limits(9) = Ide
      Itime_Limits(10) = Ihe
      Itime_Limits(11) = Ime
      Itime_Limits(12) = Ise

!  Read the first record to get some of the parameters like Lat,Long

      Call Readuf_File (Azm, Elev, Vel, Ref, MaxG, N_Gates, &
           Iy, In, Id, Ih, Im, Is, Rdel, Gtln, Rlat, Rlon, Ialt, Isn, Irn, Rname, EOF)

      If (EOF) Then
         Write (6,'("Abnormal EOF on UF file . . . ABORTING!")')
         Stop 999
      End If

      Rewind 50

      Write (LM,'(/"First UF record:"/" Azm:",F6.1/ &
          " Elev:",F6.1/                          &
          " # gates:",i6/                         &
          " Date/Time:",I4.4,"/",i2.2,"/",i2.2,1x,3i2.2/           &
          " Rdel:",F7.3/                          &
          " Gtln:",F7.3/                          &
          " Rlat:",F7.3/                          &
          " Rlon:",F8.3/                          &
          " Ialt:",i6/                            &
          " Name:",a)') Azm, Elev, N_Gates, Iy, In, Id, Ih, Im, Is, Rdel, Gtln, Rlat, Rlon, Ialt, Rname

      Isday = Id
      Xsiz = Sx * Float(Imax)
      Ysiz = Sy * Float(Jmax)
      Itime_Limits(1) = Iy
      Itime_Limits(2) = In
      Itime_Limits(3) = Id
      Itime_Limits(4) = Ih
      Itime_Limits(5) = Im
      Itime_Limits(6) = Is

!   Write the header data
  
      Dum = -999.0
      Xsize = Float(Imax) * Sx
      Ysize = Float(Jmax) * Sy
      Alt_Flag = 1.0
      Zlat = -999.0
      Zlon = -999.0
      Type = 'gnd'

      Call FDATE(Time_Current)

      Write (95) Time_Current, Imax, Jmax, Kmax, Sx, Sy, Sz, Rname, Olat, Olon, Z0, Itime_Limits, Nmosm, Init_Time, Su, Sv,  &
           Project, Dum, Dum, Dum, Namdf(1:50), Xniq, Zlat, Zlon, Type, Dum, Alt_Flag, Alt, Dum, Dum 

      Write (LM,12) Namof(1:ncho), Time_Current, Imax, Jmax, Kmax, Sx, Sy, Sz, Rname, Olat, Olon, Z0, Itime_Limits, Nmosm, &
           Init_Time, Su, Sv, Project, Namdf(1:nchd), Xniq, Rlat, Rlon, Alt, Hinf, Vint, Vslp, Xsize, Ysize, Alt_Flag, Itype 

  12  Format (/'Header for ',a,/                                    &
           '  File created on: ',a/                                 &
           '  Imax Jmax Kmax: ',3i4/                                &
           '  Sx Sy Sz [km]: ',3F6.2/                               &
           '  Radar: ',a/                                           &
           '  Olat Olon Z0 :',2F9.3,F5.1/                           &
           '  Start, End Times: ',i4.4,"/",i2.2,"/",i2.2,1x,3i2.2,1x,i4.4,"/",i2.2,"/",i2.2,1x,3i2.2/ &
           '  Nmosm:',I1,1x,'Init Time: ',i4.4,"/",i2.2,"/",i2.2,1x,3i2.2/           &
           '  Storm u,v components: ', 2f7.2/                       &
           '  Project: ',a/                                         &
           '  Name of UF input file: ',a/                           &
           '  Nyquist Interval [m/s]: ',f7.2/                       &
           '  Radar Position: ',2f8.2,' Alt [m]:',F7.1/             &
           '  Cressman Horizontal influence radius:',F5.1/          &
           '  Cressman Vertical influence:',F5.1/                   &
           '  Cressman Vertical Slope:',F5.2/                       &
           '  Domain Size:',1x,2F6.1/                               &
           '  Alt Flag:',1x,F4.1,2x,'Radar ID:',a)
  
      Call Main_Body

      Close (95)

      Call Exit
  
 1010 Write (6,'("Error:",i5," on file:",a)') Ierr, Ifile
      Call Exit

      End

      Subroutine Main_Body

      Common /One/ Imax, Jmax, Kmax
      Common /Write/ LM

      Real, Dimension (Imax,Jmax,Jmax) :: Vel, Ref, Azm, Ele
  
      Integer(Kind=2), Dimension(Imax) :: Ihd1, Ihd2, Ihd3, Ihd4
      Integer(Kind=2) Idum

      Character(Len=24) :: Time_Current

!  Initialize the arrays that will hold the cartesian data

!  Vel is radial velocity [m/s]
!  Ref is radar reflectivity [dbZ]
!  Azm is azimuth angle from North [degrees]
!  Ele is elevation angle from horizontal [degrees]

      Idum = 0

      Do i = 1, Imax
         Do j = 1, Jmax
            Do k = 1, Kmax
               Vel(i,j,k) = -999.0
               Ref(i,j,k) = -999.0
               Azm(i,j,k) = -999.0
               Ele(i,j,k) = -999.0
            End Do
         End Do
      End Do
   
!  Do the 3-D interpolation
      
      Call FDATE(Time_Current)
      Write (LM,'(/" Starting interpolation ",a/)') Time_Current

      Call Cint (Vel, Ref, Azm, Ele)

      Call FDATE(Time_Current)
      Write (LM,'(/" Stopping interpolation ",a)') Time_Current
      Call CPU_Time(CPU_secs)
      Write (LM,'(/"CPU Time used:",F10.5," seconds"/)') CPU_secs 

!  Write the interpolated data out to the disk. Use Integer*2 to save space

      Do k = 1, Kmax
         Write (LM,'(" Writing Level",i3)') k

         Do j = 1, Jmax
            Do i = 1, Imax
               Ihd1(i) = Nint(Azm(i,j,k) * 10.0)  ! =Azm from north (deg)
               Ihd2(i) = Nint(Ele(i,j,k) * 10.0)  ! =radial velocity (m/s)
               Ihd3(i) = Nint(Vel(i,j,k) * 10.0)  ! =Elevation angle (deg)
               Ihd4(i) = Nint(Ref(i,j,k) * 10.0)  ! =reflectivity (dbZ)
            End Do

            Write (95) (Ihd1(i), Ihd2(i), Ihd3(i), Ihd4(i), Idum, Idum, i=1,Imax) 
         End Do
      End Do
  
      Return
      End
 
      Subroutine Cint (VR, DZ, AZ, EL)

!     Routine to interpolate the radar spherical coordinate data to
!     Cartesian grids

      Common /One/ Imax, Jmax, Kmax
      Common /Two/ Sx, Sy, Sz, Z0, Olat, Olon, Nmosm, Su, Sv, Xniq, Alt, DZ_Name, Vr_Name
      Common /Five/ Itime_Limits(12), Init_Time(6), Isday
      Common /Write/ LM
      Common /Interp/ Hinf, Vint, Vslp
 
      Integer, Parameter :: Max_Gates = 2048

      Real, Dimension(Max_Gates) :: Ref, Vel
      Real, Dimension(Imax,Jmax,Kmax) :: WeightV, WeightZ, VR, DZ, AZ, EL

      Character(Len=8) :: Rname
      Character(Len=2) :: DZ_Name, Vr_Name

      Real :: Rah=5.885E-5, Ctr=0.0174532925
 
!     The arrays contains the grid point data.  Subroutine CINT
!     fills the grid with weighted values according to Cressman's scheme
!     with influence radius defined in terms of angles

!     The horizontal influence angle is given by Hinf and as a function
!     of range the horizontal radius of influence is 
!     R_infh = Tan(Hinf*Ctr) * Range

!     The vertical influence radius is computed as:
!         Vinf = Vint + Vslp*Elev
!         R_infv = Tan(Vinf*Ctr) * Range
!     where Vint is the vertical influence angle at Elev=0.0 deg. Vslp
!     is the change in vertical influence angle with elevation. Vslp
!     accounts for the poor vertical resolution of the VCP-12 scans
!     employed by the NWS in most storm situations. That is, the
!     influence radius grows with elevation angle.

!   The data is partitioned as follows:

!     WeightV(i,j,k), WeightZ(i,j,k) = sum of weights; Z and V are
!     handled separately because some 88D low level scans have V or Z
!     missing in alternate scans.

!      AZ (i,j,k) = weighted azimuth angle
!      EL (i,j,k) = weighted elevation angle
!      VR (i,j,k) = weighted radial velocity
!      DZ (i,j,k) = weighted reflectivity

      Logical :: TimeOnBefore, EOF

!-----------------------------------------------------------  
! Initialize the weighting fuction array

      Do i = 1, Imax
         Do j = 1, Jmax
            Do k =1, Kmax
               WeightV(i,j,k) = 0.0
               WeightZ(i,j,k) = 0.0
            End Do
         End Do
      End Do

!  Define the start, end times

      Iys = Itime_Limits(1)
      Ins = Itime_Limits(2)
      Ids = Itime_Limits(3)
      Ihs = Itime_Limits(4)
      Ims = Itime_Limits(5)
      Iss = Itime_Limits(6)
      
      Iye = Itime_Limits(7)
      Ine = Itime_Limits(8)
      Ide = Itime_Limits(9)
      Ihe = Itime_Limits(10)
      Ime = Itime_Limits(11)
      Ise = Itime_Limits(12)

! Initial time for the grid motion
  
      Iyi = Init_Time(1)
      Ini = Init_Time(2)
      Idi = Init_Time(3)
      Ihi = Init_Time(4)
      Imi = Init_Time(5)
      Isi = Init_Time(6)

      T_Init = Timz(Ihi,Imi,Isi)
      If (Idi .gt. Isday) T_Init = T_Init + 86400.0
      Stime = Timz(Ihs,Ims,Iss)
      If (Ids .gt. Isday) Stime = Stime + 86400.0
      Isn_old = -1

! Get a UF record (ray)

    1 Call Readuf_File (Azm, Elev, Vel, Ref, Max_Gates, N_Gates, Iy, In, Id, Ih, Im, Is, &
           Rdel, Gtln, Rlat, Rlon, Ialt, Isn, Irn, Rname, EOF)

      If (EOF) Then
         Write (LM,'("End of file on UF file")')
         Return
      End IF

!  Parameters read from UF file:
!     Azm - Azimuth from North
!     Elev- Elevation angle from horizontal
!     Vel - Radial velocity in each range gate
!     Ref - Radar Refectivity in each range gate
!     N_Gates - Number of range gates
!     Iy,In,Id,Ih,Im,Is - time in yymmdd hhmmss
!     Rdel - range delay (km) out to first gate
!     Gtln - gate spacing (km)
!     Rlat,Rlon - radar location
!     Ialt - radar altitude (above MSL) in m
!     Isn - sweep number
!     EOF - end of file indicator

!  Within the time window?

      If (TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iys,Ins,Ids,Ihs,Ims,Iss)) Go To 1

      If (.not.TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iye,Ine,Ide,Ihe,Ime,Ise)) Then
         Write (LM,'("End time found:",1x,i2.2,"/",i2.2,"/",i4,1x,3i2.2)') In, Id, Iy, Ih, Im, Is
         Return
      End If

      T_Init = Timz(Ihi,Imi,Isi)

!  Within the time window?

      If (     TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iys,Ins,Ids,Ihs,Ims,Iss)) Go To 1

      If (.not.TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iye,Ine,Ide,Ihe,Ime,Ise)) Then
         Write (LM,'("End time found:",1x,i2.2,"/",i2.2,"/",i4,1x,3i2.2)') In, Id, Iy, Ih, Im, Is
         Return
      End If

      If (Idi .gt. Isday) T_Init = T_Init + 86400.0
      Time  = Timz(Ih,Im,Is)
      If (Ids .gt. Isday) Time = Time + 86400.0

      Olat1 = Olat
      Olon1 = Olon

! Move the grid if needed
 
      If (Nmosm .eq. 1) Then
         Olat1 = Olat + Sv*(Time-T_Init)/ 111190.0
         Olon1 = Olon + Su*(Time-T_Init)/(111190.0*Cos(Olat*ctr))
      End If
 
!  Position of radar [km] relative to the origin

      Rx_Pos = (Rlon - Olon1) * 111.19 * Cos(Rlat*Ctr)
      Ry_Pos = (Rlat - Olat1) * 111.19

      Call Ctme(Time,Ih,Im,Is)

!  Print out some info at the beginning of each sweep

      If (Isn .ne. Isn_old) Then
         Write (LM,10) In,Id,Iy,Ih, Im, Is, Azm, Elev, Rx_Pos, Ry_Pos, Olat1, Olon1, Isn, Irn, Alt, Ialt, Z0, Rdel, &
              Gtln, Xniq, N_Gates, Rname, Time, T_Init
 
 10      Format (i2.2,"/",i2.2,"/",I4.4,1x,3i2.2,' Az:',f5.1,' El:',f5.2,&
             1x,'Xz,Yz:',2f7.1,' Corner Pt:',2f9.3,' Sweep #:',i3,&
             1x,'Rec #:',i6,&
             1x,"Alt:",F7.1,i8,' Z0:',F6.2," Rdel:",F7.3," Gtln:",F7.3,&
             1x,' Nyq:',F7.2,' Gates:',i5,1x,a,1x,"Time:",F7.0," T_init:",F7.0)
         Isn_old = Isn
      End If

!  Need twice the Nyquist interval to locally unfold data and avoid
!  averaging over aliased intervals.

      TwiceNyq = 2.0 * Xniq

!  The big loop

      Do Igate = 1, N_Gates
         If (Ref(Igate) .eq. -999.0 .and. Vel(Igate) .eq. -999.0) Cycle

!   Compute the height (z coordinate, MSL [km]) of the bin
!   and the vertical extent of the beam. Because the altitude from the
!   UF file is sometimes strange (especially for the SMART-R), use the
!   altitude from the .prm file.
 
         Range = Rdel + Float(Igate-1)*Gtln

         If (Range .le. 0.0) Cycle   ! in case the Rdel is negative

         R_infh = Tan(Hinf*Ctr) * Range
         Vinf = Vint + Vslp*Elev
         R_infv = Tan(Vinf*Ctr) * Range
         Term1 = Range * Sin(Elev*Ctr)
         Term2 = Range * Range * Rah * Cos(Elev*Ctr) * Cos(Elev*Ctr)
         z = Term1 + Term2 + Alt/1000.0
 
!   Compute the x coordinate (relative to the corner) [km]
 
         x = Rx_pos + Sin(Azm*Ctr) * Cos(Elev*Ctr) * Range
 
!   Calculate the y coordinate (relative to the corner) [km]
 
         y = Ry_pos + Cos(Azm*Ctr) * Cos(Elev*Ctr) * Range

!     To cut down on the amount of arithmetic, find the i grid points
!     within the horizontal influence radius 

         X_min = x - R_infh/2.0
         X_max = x + R_infh/2.0
         I_min = Nint(X_min/Sx) + 1
         I_max = Nint(X_max/Sx) + 1

         If (I_min .lt. 1) I_min = 1
         If (I_max .gt. Imax) I_max = Imax
         If (I_min .gt. Imax .or. I_max .lt. 1) Cycle

!     To cut down on the amount of arithmetic, find the j grid points
!     within the horizontal influence radius 

         Y_min = y - R_infh/2.0
         Y_max = y + R_infh/2.0
         J_min = Nint(Y_min/Sy) + 1
         J_max = Nint(Y_max/Sy) + 1

         If (J_min .lt. 1) J_min = 1
         If (J_max .gt. Jmax) J_max = Jmax
         If (J_min .gt. Jmax .or. J_max .lt. 1) Cycle
 
!   Find the k grid points within the influence radius
 
         Z_min = z - R_infv/2.0
         Z_max = z + R_infv/2.0
         K_min = Nint((Z_min-Z0)/Sz) + 1
         K_max = Nint((Z_max-Z0)/Sz) + 1

         If (K_min .lt. 1) K_min = 1
         If (K_max .gt. Kmax) K_max = Kmax
         If (K_min .gt. Kmax .or. K_max .lt. 1) Cycle

!   Determine weighting function for all points within the infl rad
 
         Do i = I_min, I_max
            X_pos = Float(i-1) * Sx
            Del_x = x - X_pos

            Do j = J_min, J_max
               Y_pos = Float(j-1) * Sy
               Del_y = y - Y_pos

               Do k = K_min, K_max
                  Z_pos = Float(k-1) * Sz + z0
                  Del_z = z - Z_pos
                  D_sqr = Del_x*Del_x + Del_y*Del_y + Del_z*Del_z
                  R_sqr = R_infh * R_infv

                  If (R_sqr .gt. D_sqr) Then
                     W = (R_sqr - D_sqr) / (R_sqr + D_sqr)

!  Reflectivity -  Average Z (Zed) instead of dBZ.

                     If (Ref(Igate) .ne. -999.0) Then
                        Zed = 10.0**(DZ(i,j,k) / 10.0)
                        Zed = Zed * WeightZ(i,j,k) + 10.0**(Ref(Igate) / 10.0) * W

                        If (Zed .gt. 0.0) Then
                           DZ(i,j,k) = 10.0 * ALog10(Zed/(WeightZ(i,j,k) + W))
                           WeightZ(i,j,k) = WeightZ(i,j,k) + W 
                        End If
                     End If

!     Radial Velocity - - - - - - -
!     Skip it if there is no radial velocity at the given gate

                     If (Vel(Igate) .eq. -999.0) Cycle 
                     V    = VR(i,j,k)
                     Vnew = Vel(Igate)

!     Do local "unfold" to keep from averaging over Nyquist intervals if
!     the data have not been unfolded
                     
                     If (V .ne. -999.0) Then
                        Vdiff = (V - Vnew)/TwiceNyq
                        ndif = Nint(Vdiff)
                        Vnew = Vnew + Float(ndif)*TwiceNyq   ! unfold it 
                     End If

                     If (Vnew .ne. -999.0) Then
                        V = V * WeightV(i,j,k) + Vnew * W
                        V = V / (WeightV(i,j,k) + W)
                        VR(i,j,k) = V
                     End If
!   Azimuth angle - - - - - - -
 
                     If (WeightV(i,j,k) .eq. 0.0) Then
                        S_x = 0.0
                        S_y = 0.0
                        A_sum = 0.0
                     Else
                        A_sum = AZ(i,j,k)
                        S_x = Cos(A_sum*Ctr)
                        S_y = Sin(A_sum*Ctr)
                     End If

                     A_x = Cos(Azm*Ctr)
                     A_y = Sin(Azm*Ctr)
                     Tx = S_x * WeightV(i,j,k)  + A_x * W
                     Tx = Tx / (WeightV(i,j,k) + W)
                     Ty = S_y * WeightV(i,j,k)  + A_y * W
                     Ty = Ty / (WeightV(i,j,k) + W)
                     A_sum = (Atn2(Ty,Tx)/Ctr)
                     AZ(i,j,k) = A_sum

!   Elevation angle - - - - - - -

                     A_sum = EL(i,j,k)
                     A_sum = A_sum * WeightV(i,j,k) + Elev * W
                     A_sum = A_sum  / (WeightV(i,j,k) + W)
                     EL(i,j,k) = A_sum

!   Weighting function - - - - - - -

                     WeightV(i,j,k) = WeightV(i,j,k) + W
                  End If
               End Do
            End Do
         End Do
      End Do

      Go To 1

      End
 
      Subroutine Readuf_File (Azm, Elev, VR, DZ, Max_Gates, N_Gates, Iy, In, Id, Ih, Im, Is, &
           Rdel, Gtln, Rlat, Rlon, Ialt, Isn, Irn, Rname, EOF)

!  Subroutine to read UF format tapes are return VE and DZ
!  Variables of the subroutine:

!      Man_Header  - Mandantory Header  Length = 45 words
!      Iopt_Header - Optional Header    Length = Lenoh (max of 20)
!      Loc_Header  - Local Header       Length = Lenlh (max of 20)
!      Ifld_Header - Field Header       Length = Lenfh (max of 30)
!      Idat_Header - Data Header        Length = Lendh (max of 13)
  
!  Five fields are currently supported for the SMART-R data . . . . . . .
!         Name                Field
!         ----                ------
!          VE or VR     Doppler radial velocity [m/s]
!          DM           Returned power [dBm]
!          SW           Spectral width [m/s]
!          DZ           Reflectivity [dBZ] with ground clutter cancellation
!          ZT           Raw Reflectivity [dBZ]

! For KOUN ZT is CZ
! Several other KOUN fields:

!     DR = Zdr: ratio of H and V power
!     KD = Kdp: derivative of phase
!     PH = Phase (PhiDP)
!     RH = RhoHV: correlation coefficient

      Common /Write/ LM
      Common /Two/ Sx, Sy, Sz, Z0, Olat, Olon, Nmosm, Su, Sv, Xniq, Alt, DZ_Name, Vr_Name
  
      Integer, Parameter :: Lrec = 10000         ! Nax record size in 16 bit words
      Integer, Parameter :: MaxF = 10            ! Number of fields allowed
      Integer, Parameter :: Lenoh_max = 20, Lenlh_max = 20, Lendh_max = 30
      Integer, Parameter :: Lenfh_max = 50

      Character(Len=8) :: Rname, SitNam, GenNam

      Logical :: EOF
  
      Real, Dimension(Max_Gates) :: DZ, VR

      Integer*2 Man_Header(45), &
          Iopt_Header(Lenoh_max),&
          Loc_Header(Lenlh_max),&
          Ifld_Header(Lenfh_max,MaxF),&
          Idat(Max_Gates,MaxF), Lenfh(MaxF),&
          Idat_Header(Lendh_max), Iobuf(Lrec)

      Character(Len=2) Ifname(MaxF), Char2
      Character(Len=2) :: DZ_Name, Vr_Name
      Logical :: GotDZ, GotVR

! ----------------------------------------------------------------

      EOF = .false.
      GotDZ = .false.
      GotVR = .false.

      Read (50,End=2000,Err=3000,Iostat=Ierr) Char2, Iobuf(2), (Iobuf(k),k=3,Iobuf(2)) 

!  Is this a Universal Format Tape?
  
!      Write (Char2,'(A2)') Iobuf(1)

      If (Char2 .ne. 'UF') Then
         Write (6,'("Tape not Universal Format",1x,a)') Char2
         Stop 1
      End If

!.....CHECK THAT SPECIFIED RECORD LENGTH IS LESS THAN DIM(Iobuf)=Lrec
  
      If (Iobuf(2) .gt. Lrec) Then
         Write (6,12) Iobuf(2), Lrec
  12     Format ('Record length=',I9,' exceeds buffer size=',I5)
         Stop 2
      End If
  
!.....Assign mandatory header (45 words)
  
      Do i = 1, 45
         Man_Header(i) = Iobuf(i)
      End do
 
!  The Mandatory header consists of the following format:
!   (Bull. Amer. Meteor. Soc. (1980), Vol. 61, pg. 1403)

!    Word                      Description
!    ----                      -----------
!     1                "UF" (2 char Ascii)
!     2                Record length
!     3                Position of first word of nonmandatory header
!     4                Position of first word of local header
!     5                Position of first word of data header
!     6                Record number
!     7                Volume scan number
!     8                Ray number
!     9                Record number within the ray
!    10                Sweep number
!    11-14             Radar name (ASCII)
!    15-18             Site name  (ASCII)
!    19                Degrees of latitude
!    20                Minutes of latitude
!    21                Seconds of latitude (x64)
!    22                Degrees of longitude
!    23                Minutes of longitude
!    24                Seconds of longitude (x64)
!    25                Height of antenna above sea level (meters)
!    26                Year
!    27                Month
!    28                Day
!    29                Hour
!    30                Minute
!    31                Second
!    32                Time zone (2 ASCII - "UT", "CS", "MS", etc)
!    33                Azimuth (degrees x 64)
!    34                Elevation (degrees x 64)
!    35                Sweep mode:0 - Calibration
!                                 1 - PPI
!                                 2 - Coplane
!                                 3 - RHI
!                                 4 - Vertical
!                                 5 - Target
!                                 6 - Manual
!                                 7 - Idle
!   36                 Fixed Angle (deg x 64)
!   37                 Sweep rate (deg/s x 64)
!   38                 Generation date of tape - Year
!   39                 Generation date of tape - Month
!   40                 Generation date of tape - Day
!   41-44              Generator facility name (ASCII)
!   45                 Missing or bad data flag
  
      Azm = Man_Header(33)/64.0                ! Azimuth angle
      Elev = Man_Header(34)/64.0               ! Elevation angle
      Iy = Man_Header(26)                      ! Year
      In = Man_Header(27)                      ! Month
      Id = Man_Header(28)                      ! Day
      Ih = Man_Header(29)                      ! Hour
      Im = Man_Header(30)                      ! Minute
      Is = Man_Header(31)                      ! Second

      Ibad = Man_Header(45)                    ! Bad data flag
      Isr = Man_Header(37)                     ! Sweep rate
      Igy = Man_Header(38)                     ! Year file generated
      Igm = Man_Header(39)                     ! Month file generated
      Igd = Man_Header(40)                     ! Day file generated

      Isn = Man_Header(10)                     ! Sweep number
      Irn = Man_Header(6)                      ! Record Number
      Ivn = Man_Header(7)                      ! Volume Scan Number
      Icn = Man_Header(8)                      ! Ray Number

!  Radar Name . . . . . . . . . . . . .
!  First switch the byte order of the character words

      Do i = 11, 14
         Iobuf(i) = IShftC(Iobuf(i), 8, 16)
      End Do
      
      Write (Rname,'(4a2)') (Iobuf(i),i=11,14)

!  Site Name . . . . . . . . . . . . .

      Do i = 15, 18
         Iobuf(i) = IShftC(Iobuf(i), 8, 16)
      End Do

      Write (SitNam,'(4a2)') (Iobuf(i),i=15,18)

!  Generating Facility Name . . . . . 

      Do i = 41, 44
         Iobuf(i) = IShftC(Iobuf(i), 8, 16)
      End Do

      Write (GenNam,'(4a2)') (Iobuf(i),i=41,44)

!  Correct for Y2K problems with NCAR translators

      Iy2k = 2000

      If (GenNam(1:4) .eq. 'IRIS' .or. GenNam(1:4) .eq. 'CONV' .or. GenNam(1:4) .eq. 'KOUN') Iy2k = 0
      If (GenNam(1:4) .eq. 'NCAR') Iy2k = 1900

      Iy  = Iy + Iy2k
      Igy = Igy + Iy2k    

!  Radar Latitude . . . . . . . . . . . . .

      Ideg4 = Man_Header(19)
      Imin4 = Man_Header(20)
      Isec4 = Man_Header(21)
      Rlat = Float(Ideg4) + Float(Imin4)/60.0 + Float(Isec4)/230400.0
  
!  Radar Longitude . . . . . . . . . . . . . . .

      Ideg4 = Man_Header(22)
      Imin4 = Man_Header(23)
      Isec4 = Man_Header(24)  
      Rlon = Float(Ideg4) + Float(Imin4)/60.0 + Float(Isec4)/230400.0
      Rlon = -Abs(Rlon)
  
!  Radar altitude above sea level

      Ialt = Man_Header(25)

!      Write (Char2,'(a2)') Man_Header(32)

!.....Assign Optional Header
 
      Lenoh = Man_Header(4) - Man_Header(3)

      If (Lenoh .gt. Lenoh_max) Then
         Write (6,'("Lengh of optional header too long:",i5)') Lenoh
         Stop 123
      End If

      If (Lenoh .gt. 0) Then
         Do i=1,Lenoh
            Iopt_Header(i) = Iobuf(45+i)
         End Do
      End If
  
!.....Assign Local Use header info
 
      Jstart = Man_Header(4) - 1
      Lenlh = Man_Header(5) - Man_Header(4)
 
      If (Lenlh .gt. Lenlh_max) Then
         Write (6,'("Lengh of local header too long:",i5)') lenlh
         Stop 124
      End If

      If (Lenlh .gt. 0) Then
         Do i=1, Lenlh
            Loc_Header(i) = Iobuf(Jstart+i)
         End Do
      End If
 
!.....Assign data header info
 
      Jstart = Man_Header(5)
      Nfe = Iobuf(Jstart)       ! Number of data fields

      If (Nfe .gt. MaxF) Then
         Write (6,'(i2," data fields currently allowed - need ",i5)') MaxF, Nfe 
         Stop 3
      End If
 
      Lendh = 3 + 2*Nfe

      If (Lendh .gt. Lendh_max) Then
         Write (6,'("Length of data header too long:",i5)') Lendh
         Stop 126
      End If

      Do j=1, Lendh
         Idat_Header(j) = Iobuf(Jstart + j - 1)
      End Do
 
!.....ASSIGN FIELD NAME AND FIELD HEADER INFO

!  Field Header Format is as follows:
!     Word                       Description
!     ----                       ------------
!      1                  Position of first data word
!      2                  Scale factor
!      3                  Range to first gate [km]
!      4                  Adjustment to center of first gate [m]
!      5                  Sample volume spacing [m]
!      6                  Number of sample volumes
!      7                  Sample volume depth [m]
!      8                  Horizontal beam width [deg x 64]
!      9                  Vertical beam width [deg x 64]
!     10                  Reciever bandwidth [Mhz]
!     11                  Transmitted Polarization: 0 = horizontal
!                                                   1 = vertical
!                                                   2 = circular
!                                                  >2 = elliptical
!     12                  Wavelength [cm x 64]
!     13                  Number of samples used in field estimate
!     14                  Threshold field [2 ASCII]
!     15                  Threshold value
!     16                  Scale
!     17                  Edit code [2 ASCII]
!     18                  PRT [microseconds]
!     19                  Bits per sample volume [generally 16]
!     20-?                Individual field words, as follows:
!           For VE
!           ------
!     20                  Nyquist velocity [scaled]
!     21                  "FL" [2 ASCII] if flagged in least sig bit
!                          with NCAR bad velocity flag [1=good,0=bad]
!           For DM
!           ------
!     20                  Radar Constant
!     21                  Noise power  [dBmW x scale]
!     22                  Reciever gain [dB x scale]
!     23                  Peak Power [dBmW x scale]
!     24                  Antenna Gain [dB x scale]
!     25                  Pulse duration [microsec x 64]

      Do j = 1, Nfe
 
!        GET FIELD NAME, but first switch the byte order

         Idat_Header(2+(2*J)) = IShftC(Idat_Header(2+(2*J)), 8, 16) 
         Write (Ifname(j),'(a2)') Idat_Header(2+(2*J))

!        GET STARTING WORD AND LENGTH FOR THIS FIELD
 
         Jstart = Idat_Header(3+(2*J))
         Lenfh(J) = Iobuf(Jstart) - Jstart
         Jstart = Jstart - 1

         If (Lenfh(j) .gt. Lenfh_max) Then
            Write (6,'("Length of field header too big:",i5)') Lenfh(j)
            Stop 125
         End If
 
!   Define Field Header
  
         Do I = 1, Lenfh(j)
            Ifld_Header(i,j) = Iobuf(Jstart+I)
         End Do
      End Do

!........GET STARTING POINTER FOR DATA IN EACH FIELD
 
      Do j = 1, Nfe
         Jstart = Ifld_Header(1,j) - 1
         N_gates = Ifld_Header(6,j)

!  Load up only as many gates as the main program has allocated for

         If (N_Gates .gt. Max_Gates) Then
            N_Gates = Max_Gates
         End If

         Do i = 1, N_gates
            Idat(i,j) = Iobuf(Jstart+i)
         End Do
      End Do

      xnyq2 = -999.0
      xnyq1 = -999.0
      n_gatesr = -999
      adjr = -999.0
      gtr = -999.0
      rdelr = -999.0

      Do j = 1, Nfe

!  Find the reflectivity data

         If (Ifname(j) .eq. DZ_Name) Then
            GotDZ = .true.
            Scale = Ifld_Header(2,j)
            RdelR = Ifld_Header(3,j)
            GtR   = Ifld_Header(5,j)
            GtR   = GtR/1000.0
            AdjR  = Ifld_Header(4,j)
            AdjR  = AdjR/1000.0
            N_gatesR = Ifld_Header(6,j)
            RangeOne = RdelR + AdjR        ! range of gate one
            ij = 0
            If (N_GatesR .gt. Max_Gates) N_GatesR = Max_Gates

            Do i = 1, N_gatesR
               Range= Float(i-1) * GtR + RangeOne

               If (Range .gt. 0.0) Then     ! only use positive ranges
                  ij = ij + 1
                  Z = Idat(i,j)
                  DZ(ij) = Z/Scale
                  If (Idat(i,j) .eq. Ibad) DZ(ij) = -999.0
               End If
            End Do

! Readjust some values in case first gate(s) discarded because of (-) range

            RdelR = RdelR + Float(N_GatesR -ij) * GtR ! Readjust RdelR
            N_GatesR = ij                             ! Readjust N_gatesR
         Else If (Ifname(j) .eq. Vr_Name) Then        !  Get Velocity data  
            GotVR = .true.
            Scale = Ifld_Header(2,j)
            RdelV = Ifld_Header(3,j)
            IgtV   = Ifld_Header(5,j)
            GtV = Float(IgtV)/1000.0
            AdjV  = Ifld_Header(4,j)
            AdjV = AdjV/1000.0
            N_GatesV = Ifld_Header(6,j)
            N_GatesVOrig = N_GatesV
            Wavelength = Ifld_Header(12,j)
            Wavelength = Wavelength/6400.0
            Prt = Ifld_Header(18,j)
            Prf = (1.0/Prt)*1000000.0
            Xnyq1 = Prf*Wavelength/4.0
            Xnyq2 = Ifld_Header(20,j)
            Xnyq2 = Xnyq2/Scale
            RangeOne = RdelV + AdjV        ! range of gate one
            ij = 0
            If (N_GatesV .gt. Max_Gates) N_GatesV = Max_Gates

            Do i = 1, N_gatesV
               Range = Float(i-1)* GtV + RangeOne

               If (Range .gt. 0.0) Then
                   ij = ij+1
                   VelR = Idat(i,j)
                   VR(ij) = VelR/Scale
                   If (Idat(i,j) .eq. Ibad) VR(ij) = -999.00
               End If
            End Do               ! i=1,ngates

! Readjust some values in case first gate(s) discarded because of (-) range

            RdelV = RdelV + Float(N_GatesV -ij) * GtV ! Readjust RdelV
            N_GatesV = ij                             ! Readjust N_gatesV
         End If
      End Do                     ! j=1,Nfe

      If (GotVR .and. GotDZ) Then
         N_Gates = N_GatesR
         Gtln = GtR
         Adj = AdjR
         Rdel = RdelR + Adj
         GtlnV = GtV
         GtlnR = GtR
         RdelV = RdelV + AdjV
         RdelR = RdelR + AdjR
         Xniq = Xnyq1
         Return
      Else
         Write (6,'("One of the fields ",a,1x,a, " is missing!")') DZ_Name, VR_Name
         Write (6,'("Fields available ",10(a,1x))') (Ifname(j),j=1,NFE)
         Stop 'missing fields'
      End If

 1020 Write (6,'("Error= ",i5," on file: ",a)') Ierr, File
      Stop 1020

 2000 EOF = .True.
      Close (50)

      Return

 3000 Write (6,'("Error on UF read: ",i5)') Ierr

!     The most common error here from g77 is 110 "Reading past end of
!     record"

      Stop 3000

      End

      Block Data

      Common /One/ Imax, Jmax, Kmax
      Common /Two/ Sx, Sy, Sz, Z0, Olat, Olon, Nmosm, Su, Sv, Xniq, Alt, DZ_Name, Vr_Name
      Common /Five/ Itime_Limits(12), Init_Time(6), Isday
      Common /Write/ LM
      Common /Interp/ Hinf, Vint, Vslp

      Character(Len=2) :: DZ_Name, Vr_Name

      Data LM /60/
  
      End
