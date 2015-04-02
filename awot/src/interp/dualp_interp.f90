Program dualp_interp

  !--------------------------------------------------------------
  !     Program to interpolate ground-based dual pol radar date to 2-D grid using
  !     a "hybrid" scan produced by a DEM topo map. A hybrid scan is a
  !     polar coordinate table that details which tilt level to use at
  !     each range bin so that the beam clears the topo all the way back
  !     to the radar. Hybrid reflectivity is a horizontal map of all the
  !     hybid bins. Once the spherical coordinate data is placed on the
  !     polar grid it is interpolated to Cartesian coordinates using a
  !     Cressman interpolation scheme.
  !--------------------------------------------------------------

  ! Based on the hybrid_interp program for constructing hybrid reflectivity
  
  !  David Jorgensen NOAA/NSSL May 2010 for CPOL SMART-R
  !  All printed output is directed to a .log file
  !  Radar data is input as Universal Format (UF)
  
  !  The parameters used in the program:
  
  !  Namof- name of the disc file that will contain output Cartesian data
  !  Namdf- name of the disc file that contains the UF 88D data
  !  Imax - maximum points in X direction (East-West)
  !  Jmax - maximum points in Y direction (North-South)
  !  Sx   - distance between points (km) in X direction
  !  Sy   - distance between points (km) in Y direction
  
  !  Nmosm- flag for type of grid reference desired
  !       - =0 with respect to the ground; in this case,
  !          OLAT - latitude of lower left hand corner
  !          OLON - longitude of lower left hand corner
  !       - =1 with respect to the moving grid center; in that case,
  !          Olat = Olat + Sv*(Time-T_init)/111120.0
  !          OLon = Olon + Su*(Time-T_init)/(111120.0*Cos(Olat))
  !  IHS, IMS, ISS - starting time
  !  IHE, IME, ISE - ending time
  
  !  Summary of disc lu's used:
  !       10 - parameter file
  !       64 - log file
  !       50 - input UF file
  !       60 - input hybrid scan file
  !       66 - input topo file (for VPR correction)
  !       95 - OUTPUT (Cartesian) file
  
  !     The input files are on the command line to facilitate scripting
  !     the entire IOP
  !--------------------------------------------------------------
  
  Integer FirstChar
  
  Integer, Parameter :: MaxG = 2048

  Real, Dimension(MaxG) ::  VR, DZ, ZT, DR, KD

  Logical OKTime, EOF

  Common /One/ Imax, Jmax
  Common /Two/ Sx, Sy, Olat, Olon, Nmosm, Su, Sv, Xnyq, Alt, dbz_offset, Ivpr, Iatten
  Common /Thr/ Nazms, Nbins, Azm_Step, RdelH, GtlnH, RlatH, RlonH
  Common /For/ Ioffset, Sectors(2,25), NOffsets(25), MaxS(25), Idbz
  Common /Five/ Itime_Limits(12), Init_Time(6), Isday, Iadjyr
  Common /Write/ LM
  Common /Interp/ Hinf
  Common /CPU/ Machine
  Common /VRP/ DelZ1, tanb
  Common /Proj/ Project
  Common /Freezing/ Fr_Hgts(1000,2), Nfrz, Freezing_File
  
  Character(Len=100) Line, Ifile, PrmFile, Namdf, Namof, LogFile, HyFile, TopoF, Freezing_File
  Character Parms(4)*6
  
  Character Time_Current*24, UserN*8, Itype*3, Ident*8
  Character Rname*8, Project*16, Project1*16, Desc*57
  
  Data Parms(1) /'060101'/, Parms(2) /'000000'/, Parms(3) /'101231'/, Parms(4) /'235959'/
  Data Ident /'hybdualp'/
  
  !--------------------------------------------------------------
  !  Dynamic memory for the big arrays
  !--------------------------------------------------------------

  Real, Allocatable :: Cart(:,:,:)     ! (x,y)
  Real, Allocatable :: Hybrid(:,:)     ! (azms,bins)
  Real, Allocatable :: Polar(:,:,:)    ! (azms,bins,4) 1=dbz; 2=elev; 3=Zdr; 4=Kdp
  Real, Allocatable :: Topo(:,:)       ! (x,y)
  
  !--------------------------------------------------------------
  !    Open the data file containing the parameter information.
  !--------------------------------------------------------------
  
  nparms = IargC()
  
  If (nparms .lt. 2) Then
     Write (6,'("Enter prm file name:",$)')
     Read (5,'(a)') PrmFile
     Write (6,'("Enter UF file name:",$)')
     Read (5,'(a)') Namdf
  Else
     Call GetArg(1,PrmFile)
     Call GetArg(2,Namdf)
     
     IF (nparms .gt. 2) Then
        Do i = 1, nparms
           Call GetArg(i+2, Parms(i))
        End Do
     End If
  End If
  
  !--------------------------------------------------------------
  !     Which machine are we running on? Makes a difference in how the
  !     ASCII characters are imported. On intel macs the ASCII characters
  !     are switched
  !--------------------------------------------------------------
  
  Machine = 0
  Call GetEnv('MACHINE',Line)
  IF (Line(1:1) .eq. 'i') Machine = 1
  
  !--------------------------------------------------------------
  !     If the start & end times were on the command line, then decode
  !     them. If not then dummy up the stimes to use all the sweeps in the
  !     UF file
  !--------------------------------------------------------------
  
  Read (Parms(1),'(3i2)') Iys, Ins, Ids
  Read (Parms(2),'(3i2)') Ihs, Ims, Iss
  Read (Parms(3),'(3i2)') Iye, Ine, Ide
  Read (Parms(4),'(3i2)') Ihe, Ime, Ise
  
  nchp = LenTrim(PrmFile)
  Ifile = PrmFile(1:nchp)
  
  Open (10, File=Ifile, Err=1010, Iostat=Ierr, Status='Old')
  
  !--------------------------------------------------------------
  !     Open the UF data file
  !--------------------------------------------------------------
  
  nchd = LenTrim(Namdf)
  Ifile = Namdf(1:nchd)
  
  Open (50,File=Ifile, Form='UNFORMATTED', Err=1010, Iostat=Ierr, Status='Old', Convert='Big_Endian')
  
  !--------------------------------------------------------------
  !     Name the log file 'dplyymmddhhmmss.log Strip off the
  !     directory info so the .log file will be written in the current
  !     directory.
  !--------------------------------------------------------------
  
  nchl = FirstChar(Namdf,'.') - 1
  nchf = LastChar(Namdf,'/') + 1
  LogFile = Namdf(nchf:nchl) // '.log'
  LogFile(1:3) = 'dpl'
  
  !--------------------------------------------------------------
  !  Check to see if the file name is from the NCAR translator UF file or the SIGMET UF file format
  !  The difference is in the file names: NCAR: ufd.1051130211324.RADAR_NAME.0.tape
  !                                     SIGMET: SR1071219010005.RAWJ5M2.UF
  !--------------------------------------------------------------

  nchu = LastChar(Namdf,'.') + 1
  
  If (Namdf(nchu:) .eq. 'tape') Then
     nchu = LastChar(Namdf,'S') - 1
     Logfile(1:3) = 'dpl'             ! prefix "dpl" for dual-pol parameters
     Logfile(4:) = namdf(nchf+5:nchu)
     
     nch = LenTrim(LogFile)
     LogFile = Logfile(1:nch) // 'log'
  End If
  
  nchl = LenTrim(LogFile)
  Ifile = Logfile(1:nchl)
  
  If (LM .ne. 6) Then
     Open (LM, File=Ifile, Iostat=Ierr, Err=1010)
  End If
  
  !--------------------------------------------------------------
  !  Write out info to the .log file, including the user name
  !--------------------------------------------------------------
  
  UserN = 'xxxxxxxx'
  Call GetLog(UserN)
  
  Write (LM,'("Hybrid-Scan Radar Interpolation program - .prm file:",a, &
       " for user: ",a)') Prmfile(1:nchp), UserN
  Write (LM,'("  Log File:",a)') LogFile(1:nchl)
  
  !------------------------------------------------------------
  !     Name the output file as the .log file but append ".dat" as
  !     the extension and put it in the same directory as the .log file
  !------------------------------------------------------------
  
  nch = LenTrim(LogFile)
  Namof = Logfile(1:nch-3) // 'dat'
  
  ncho = LenTrim(Namof)
  Ifile = Namof(1:ncho)
  
  Open (95, File=Ifile, Iostat=Ierr, Err=1010, Form='Unformatted')
  Write (LM,'(/"Output file name: ",a)') Namof(1:ncho)
  
  !--------------------------------------------------------------      
  !  Read parameters from the .prm file
  !     The first line is the indicator for the type of .prm file
  !     (hybrid_interp.prm)
  !--------------------------------------------------------------
  
  Read (10,'(a)') Line
  Call CaseFold(Line)
  nch = FirstChar(Line,' ') - 1
  
  If (Line(1:nch) .ne. 'DUALP_INTERP') Then
     Write (LM,'(/"  Improper type of .prm file - must start with ", &
          "dualp_interp ",a/)') Line(1:nch)
     Call Exit
  End If
  
  Read (10,'(a)') Project
  
  !--------------------------------------------------------------
  !     Read the hybrid scan file. A hybrid scan file is in polar
  !     coordinates (theta,range) that matches the input radar data
  !     resolution. For the SMART-R in HMT it was every 1 deg and 133 m in
  !     range. For Debris Flow in 2007 the resolution was 1 deg, 100 m
  !--------------------------------------------------------------
  
  Read (10,'(a)') HyFile
  
  nch = FirstChar(HyFile,' ') - 1
  HyFile = HyFile(1:nch)
  Ifile = HyFile
  
  Open (60, File=Ifile, Iostat=Ierr, Err=1010, Form='Unformatted')
  
  Write (LM,'("Hybrid Scan File Name= ",a)') HyFile(1:nch)
  
  Read (60) Time_Current, Desc, Rname, Project1, Rlath, Rlonh, Nbins, Nazms, RdelH, GtlnH, Azm_Step
  
  Write (LM,'(/"Hybrid Scan parameters:"/ &
       " File Created ",a/ &
       " Description: ",a/ &
       " Radar: ",a/ &
       " Project: ",a/ &
       " Radar Pos:",1x,F8.3,1x,F9.3/ &
       " Num bins: ",i4/ &
       " Num azms: ",i4/ &
       " Range Delay:",F8.3/ &
       " Gate Length:",F8.3/ &
       " Azm Step:",F6.1/)') Time_Current, Desc, Rname, Project, Rlath,&
       Rlonh, Nbins, Nazms, RdelH, GtlnH, Azm_Step 
  
  !  The grid parameters
  
  Read (10,*) Imax, Jmax, Sx, Sy
  Read (10,*) Nmosm, Olat, Olon, Iynd, Ihms, Su, Sv
  
  Call Tcnvt(Iynd, Iyi, Ini, Idi)
  Call Tcnvt(Ihms, Ihi, Imi, Isi)
  
  !--------------------------------------------------------------
  !     The altitude (MSL) of the radar. Need this becasue sometimes the
  !     altitude from the UF file is weird, especially for the SMART-R.
  !     dbz_offset is the reflectivity bias
  !     Ivpr is whether to do a VPR correction (0=No, 1=Yes)
  !     If VPR is desired, then need a topo file
  !     Iatten is do the attenuation correction (1=yes 0=no)
  !--------------------------------------------------------------
  
  Read (10,*) Alt, dbz_offset, Ivpr, Iatten
  
  DelZ1 = -999.0
  tanb = -999.0
  
  If (Ivpr .eq. 1) Then
     Read (10,*) DelZ1, tanb
     Read (10,'(a)') TopoF
     Read (10,'(a)') Freezing_File
     
     !     Take out the comments after the first blank character
     
     nch = FirstChar(Freezing_File,' ') - 1
     Freezing_File = Freezing_File(1:nch)
  End If
  
  !--------------------------------------------------------------
  !  Cressman horizontal weighting interval, elev offset flag, & dbz flag
  !  Idbz = 0 uses the processed or filtered reflectivity
  !  Idbz = 1 uses the raw (ZT) field
  !  Ioffset is the number of azimuth sectors
  !  Noffsets(Ioffset) is the number of elevation steps to add in each sector
  !    to account for local beam blockages not taken care of in the hybrid scan 
  !    design
  !  MaxS is the maximum step increase
  !--------------------------------------------------------------
  
  Read (10,*) Hinf, Ioffset, Idbz
  
  If (Ioffset .gt. 25) Then
     Write (6,'("Too many azimuth sectors, must by <26")')
     Call Exit
  End If

  If (Ioffset .gt. 0) Then
     Do i = 1, Ioffset
        Read (10,*) Sectors(1,i), Sectors(2,i), NOffsets(i), MaxS(i)
        Write (LM,'(" Az sector from: ",F5.1," to ",F5.1, &
             " elev step offset: ",i2," max step:",i3)') &
             Sectors(1,i), Sectors(2,i),NOffsets(i), MaxS(i) 
     End Do
  Else
     Write (LM,'(/"No offsets added to hybrid-scan tables")')
  End If
  
  Close(10)
  
  Write (LM,'(/"Input Variables:"/ &
       " Project Name:",1x,a/ &
       " Hybrid File: ",1x,a/ &
       " Imax:        ",1x,i4/ &
       " Jmax:        ",1x,i4/ &
       " Sx:          ",1x,F6.3/ &
       " Sy:          ",1x,F6.3/ &
       " Nmosm:       ",1x,I2/ &
       " Olat:        ",1x,F7.3/ &
       " Olon:        ",1x,F8.3/ &
       " Init Time:   ",1x,i6,1x,i6/ &
       " Su:          ",1x,F6.2/ &
       " Sv:          ",1x,F6.2/ &
       " Alt:         ",1x,F6.2/ &
       " dBZ offset:  ",1x,F6.2/ &
       " VPR Corr:    ",1x,i1/   &
       " Atten Corr:  ",1x,i1/   &
       " Del_dbz1:    ",1x,F7.2/ &
       " tan B:       ",1x,F7.2/ &
       " Idbz:        ",1x,i2/ &
       " Hinf:        ",1x,F7.2/)') &
       Project, HyFile, Imax, Jmax, Sx, Sy, Nmosm, &
       Olat, Olon, Iynd, Ihms, Su, Sv, Alt, dbz_offset, Ivpr, Iatten, &
       Delz1, tanb, Idbz, Hinf
  
  If (Idbz .eq. 0) Then
     Write (LM,'("using DZ field for Reflectivity")')
  Else
     Write (LM,'("using ZT field for Reflectivity")')
  End If
  
  If (Ivpr .eq. 1) Then
     Write (LM,'(/"Freezing level file:",a)') Freezing_File
  End If
  
  Iadjyr = 1900
  If (Iyi .lt. 75) Iadjyr = 2000
  
  !  Check the grid motion init time
  
  If (Nmosm .eq. 1) Then
     If (.not. OKTime(Iyi+Iadjyr,Ini,Idi,Ihi,Imi,Isi)) Then
        Write (LM,'("Time: ",i2.2,"/"i2.2,"/",i2.2,1x,3i2.2 ," is not OK")') Iyi, Ini, Idi, Ihi, Imi, Isi
        Call Exit
     End If
  End If
  
  Init_Time(1) = Iyi
  Init_Time(2) = Ini
  Init_Time(3) = Idi
  Init_Time(4) = Ihi
  Init_Time(5) = Imi
  Init_Time(6) = Isi
  
  !  Are the start times OK?
  
  If (.not. OKTime(Iys+Iadjyr,Ins,Ids,Ihs,Ims,Iss)) Then
     Write (LM,'("Time: ",i2.2,"/"i2.2,"/",i2.2,1x,3i2.2," is not OK")') Iys,Ins,Ids,Ihs,Ims,Iss
     Stop
  End If
  
  !  Are the end times OK?
  
  If (.not. OKTime(Iye+Iadjyr,Ine,Ide,Ihe,Ime,Ise)) Then
     Write (LM,'("Time: ",i2.2,"/"i2.2,"/",i2.2,1x,3i2.2," is not OK")') Iye,Ine,Ide,Ihe,Ime,Ise
     Call Exit
  End If
  
  Itime_Limits(1) = Iys
  Itime_Limits(2) = Ins
  Itime_Limits(3) = Ids
  Itime_Limits(4) = Ihs
  Itime_Limits(5) = Ims
  Itime_Limits(6) = Iss
  Itime_Limits(7) = Iye
  Itime_Limits(8) = Ine
  Itime_Limits(9) = Ide
  Itime_Limits(10) = Ihe
  Itime_Limits(11) = Ime
  Itime_Limits(12) = Ise
  
  !  Read the first record to get some of the parameters like Lat,Long
  
  Call Readuf_File (Azm, Elev, VR, DZ, ZT, DR, KD, MaxG, N_Gates, Iy, In, Id, Ih, Im, Is, Rdel, Gtln, Rlat, Rlon, Ialt, & 
       Isn, Irn, Rname, EOF)
  
  If (EOF) Then
     Write (6,'("Abnormal EOF on UF file . . . ABORTING!")')
     Call Exit
  End If
  
  Rewind 50
  
  Write (LM,'(/"First UF record:"/" Azm:",F6.1/ &
       " Elev:",F6.1/ &
       " # gates:",i9/ &
       " Date/Time:",3i2.2,1x,3i2.2/ &
       " Rdel:",F7.3/ &
       " Gtln:",F7.3/ &
       " Rlat:",F12.8/ &
       " Rlon:",F13.8/ &
       " Ialt:",i6/ &
       " Name:",a)') Azm, Elev, N_Gates, Iy,In,Id,Ih,Im,Is, Rdel, Gtln, Rlat, Rlon, Ialt, Rname
  
  !  Make sure the radars match
  
  If (Abs(Rlat-RlatH) .gt. 0.001 .or. Abs(Rlon-RlonH) .gt. 0.001) Then 
     Write(LM,'("Radar Locatioms do not match . . aborting!")')
     Write (LM,*) 'Lats: ', Rlat, RlatH, ' Lons: ', Rlon, RlonH
     Call Exit
  End If
  
  If (Rdel .ne. RdelH .or. Gtln .ne. GtlnH) Then
     Write(LM,'("Radar parms do not match . . aborting!")')
     Write (LM,*) 'Rdel: ', Rdel,RdelH, ' Gtln: ', Gtln, GtlnH
     Call Exit
  End If
  
  Isday = Id
  Itime_Limits(1) = Iy
  Itime_Limits(2) = In
  Itime_Limits(3) = Id
  Itime_Limits(4) = Ih
  Itime_Limits(5) = Im
  Itime_Limits(6) = Is
  
  !   Write the header data
  
  Dum = -999.0
  Itype = 'gnd'
  Xsize = Float(Imax) * Sx
  Ysize = Float(Jmax) * Sy
  Alt_Flag = 1.0
  
  Call FDATE(Time_Current)
  
  Write (95) Ident, Time_Current, Imax, Jmax, Sx, Sy, Rname, &
       Olat, Olon, Itime_Limits, Nmosm, Init_Time, &
       Su, Sv, Project, Ivpr, DelZ1, tanb, Namdf, &
       Xnyq, Rlat, Rlon, Itype, Alt_Flag, Alt, dbz_offset, Iatten, Dum
  
  Write (LM,12) Namof(1:ncho), Time_Current, Imax, Jmax, Sx, &
       Sy, Rname, Olat, Olon, Itime_Limits, Nmosm, &
       Init_Time, Su, Sv, Project, Namdf(1:nchd), Xnyq, Rlat, Rlon, &
       Alt, dbz_offset, Hinf, Xsize, Ysize, Alt_Flag, Ivpr, Iatten, DelZ1, tanb 
  
12 Format (/'Header for ',a,/ &
        '  File created on: ',a/ &
        '  Imax Jmax: ',2i4/ &
        '  Sx Sy [km]: ',2F6.2/ &
        '  Radar: ',a/ &
        '  Olat Olon:',2F9.3/ &
        '  Start, End Times: ',3i2.2,1x,3i2.2,1x,3i2.2,1x,3i2.2/ &
        '  Nmosm:',I1,1x,'Init Time: ',3i2.2,1x,3i2.2/ &
        '  Storm u,v components: ', 2f7.2/ &
        '  Project: ',a/ &
        '  Name of UF input file: ',a/ &
        '  Nyquist Interval [m/s]: ',f7.2/ &
        '  Radar Position: ',2f13.7,' Alt [m]:',F6.1/ &
        '  dBZ offset:',F6.2/ &
        '  Cressman Horizontal influence radius:',F5.1/ &
        '  Domain Size:',1x,2F6.1/ &
        '  Alt Flag:',1x,F4.1/ &
        '  IVPR:',i1/ &
        '  Iatten:',i1/ &
        '  DelZ1:',F7.2/ &
        '  tanb:',F7.2)
  
  Allocate (Polar(Nazms,Nbins,4), Hybrid(Nazms,Nbins))
  
  Do na = 1, Nazms
     Read (60) (Hybrid(na,nb),nb=1,Nbins)
  End Do
  
  Close (60)
  
  !  Read in the topo data
  !  If VPR is desired, read in the melting level heights
  
  If (IVPR .eq. 1) Then
     Allocate (Topo(Imax,Jmax))
     Call RdTopo(TopoF, Topo)
     Call Load_Freezing_Level
  End If
  
  !  Compute the hybrid reflectivity on a polar grid
  
  Call Hreflec (Hybrid, Polar, Topo)
  
  Deallocate(Hybrid)
  If (IVPR .eq. 1) Deallocate(Topo)
  
  Allocate (Cart(Imax,Jmax,3))
  
  !  Fill up the missing values in Polar (-998.0) with -999.0
  !  Do the interpolotion in Z space
  
  Do na = 1, Nazms
     Do nb = 1, Nbins
        If (Polar(na,nb,1) .eq. -998.0) Polar(na,nb,1) = -999.0
        If (Polar(na,nb,2) .eq. -998.0) Polar(na,nb,2) = -999.0
        If (Polar(na,nb,1) .ne. -999.0) Polar(na,nb,1) = 10.0**(Polar(na,nb,1)/10.0)
     End Do
  End Do
  
  !  Do the 3-D interpolation
  
  Call Cnvt_Cartesian (Polar, Nazms, Nbins, RdelH, GtlnH, Azm_Step, &
       Cart, Imax, Jmax, Sx, Sy, Olat, Olon, Rlath, Rlonh, Hinf) 
  
  !  Convert back to dBZ
  
  Do i = 1, Imax
     Do j = 1, Jmax
        If (Cart(i,j,1) .ne. -999.0) Cart(i,j,1) = 10.0 * Alog10(Cart(i,j,1))
     End Do
  End Do
  
  Deallocate(Polar)
  Call FDATE(Time_Current)
  Write (LM,'(" Stopping interpolation ",a)') Time_Current
  Call CPU_Time(CPU_secs)
  Write (LM,'(/"CPU Time used:",F10.5," seconds"/)') CPU_secs 
  
  !  Write the interpolated data out to the disk
  
  Write (LM,'(" Writing Output File"/)')
  
  Do j = 1, Jmax
     Write (95) (Cart(i,j,1),i=1,Imax) 
     Write (95) (Cart(i,j,2),i=1,Imax) 
     Write (95) (Cart(i,j,3),i=1,Imax) 
  End Do
  
  Deallocate(Cart)
  
  Close (95)
  
  Call Exit
  
1010 Write (6,'("Error:",i5," on file:",a)') Ierr, Ifile
  Call Exit
  
End Program dualp_interp

Subroutine Hreflec(Hybrid, Polar, Topo)
  
  !--------------------------------------------------------------
  !     Routine to interpolate the radar spherical coordinate data to
  !     a hybrid scan reflectivity using a lookup of elevation angles
  !     Special sectors are allowed in which the elevation angles are
  !     increased to get around local beam blockage by close obstacles
  !--------------------------------------------------------------
  
  Common /One/ Imax, Jmax
  Common /Two/ Sx, Sy, Olat, Olon, Nmosm, Su, Sv, Xnyq, Alt, dbz_offset, Ivpr, Iatten
  Common /Thr/ Nazms, Nbins, Azm_Step, RdelH, GtlnH, RlatH, RlonH
  Common /For/ Ioffset, Sectors(2,25), NOffsets(25), MaxS(25), Idbz
  Common /Five/ Itime_Limits(12), Init_Time(6), Isday, Iadjyr
  Common /Write/ LM
  Common /Interp/ Hinf
  
  Integer, Parameter :: Nelevs = 14, Max_Gates = 2048
  
  Real, Dimension(Max_Gates) :: Ref, Vel, ZT, DR, KD
  Real, Dimension(Nazms,Nbins,4) :: Polar
  Real, Dimension(Nazms,Nbins) :: Hybrid
  Real, Dimension(Imax,Jmax) :: Topo
  Integer, Dimension(Nelevs) :: nVCP12

!  Dimension Polar(Nazms,Nbins,2), Hybrid(Nazms,Nbins), Topo(Imax,Jmax), nVCP12(Nelevs)
  
  Character Rname*8
  Logical TimeOnBefore, EOF
  
  Data Ctr /0.01745329/

  ! Tilts for VCP-12 
  
  Data nVCP12 /5, 9, 13, 18, 24, 31, 40, 51, 64, 80, 100, 125, 156, 195/

  !--------------------------------------------------------------
  !  Initilize the output array and
  !  Convert the Hybrid scan into elev step numbers for VCP-12
  !     Polar(iazm,igate,1) contains the hybrid scan reflectivity. It is
  !     initialized to -998.0 so that missing reflectivity (-999.0) can be
  !     a valid fill and high elevation data won't be blended with low
  !     -elev data if there is no echo at low levels.
  !--------------------------------------------------------------

  Do na = 1, Nazms
     Do nb = 1, Nbins
        Polar(na,nb,1) = -998.0       
        Polar(na,nb,2) = -998.0
        nElev = Nint(Hybrid(na,nb)*10.0)
        nstep = -999
        
        !  Clean up the elevs so that they match the vcp elevs
        
        Do ne = 1, Nelevs
           If (nElev .eq. nVCP12(ne)) nstep = ne
        End Do
        
        If (nstep .gt. 0) Then
           Hybrid(na,nb) = nstep
        Else
           Hybrid(na,nb) = -999.0
        End If
     End Do
  End Do

  !  Define the start, end times
  
  Iys = Itime_Limits(1) + Iadjyr
  Ins = Itime_Limits(2)
  Ids = Itime_Limits(3)
  Ihs = Itime_Limits(4)
  Ims = Itime_Limits(5)
  Iss = Itime_Limits(6)
  
  Iye = Itime_Limits(7) + Iadjyr
  Ine = Itime_Limits(8)
  Ide = Itime_Limits(9)
  Ihe = Itime_Limits(10)
  Ime = Itime_Limits(11)
  Ise = Itime_Limits(12)
  
  !  Initial time for the grid motion
  
  Iyi = Init_Time(1) + Iadjyr
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
  Ntotal_Azms = -999
  total_vpr_corr = 0.0
  ntotal_vpr = 0
  
  !  Get a UF record (ray)
  
1 Call Readuf_File (Azm, Elev, Vel, Ref, ZT, DR, KD, Max_Gates, N_Gates, &
       Iy, In, Id, Ih, Im, Is, Rdel, Gtln, Rlat, Rlon, Ialt, Isn, Irn, Rname, EOF)
  
  If (EOF) Then
     Write (LM,'("End of file on UF file")')
     
     If (Ivpr .eq. 1) Then
        total_vpr_corr = total_vpr_corr/Float(ntotal_vpr)
        Write (LM,'(/"Ave vpr correction:",1x,F9.2,1x,"N=",i8)') total_vpr_corr, ntotal_vpr
     End If
     
     Polar(1,1,1) = 100.0
     
     Return
  End If
  
  !--------------------------------------------------------------
  !  Parameters read from UF file:
  !     Azm - Azimuth from North
  !     Elev- Elevation angle from horizontal
  !     Vel - Radial velocity in each range gate
  !     Ref - Radar Refectivity (filtered) in each range gate
  !     ZT  - Raw Radar Reflectivity in each gate
  !     N_Gates - Number of range gates
  !     Iy,In,Id,Ih,Im,Is - time in yymmdd hhmmss
  !     Rdel - range delay (km) out to first gate
  !     Gtln - gate spacing (km)
  !     Rlat,Rlon - radar location
  !     Ialt - radar altitude (above MSL) in m
  !     Isn - sweep number
  !     EOF - end of file indicator
  
  !  Within the time window?
  !--------------------------------------------------------------
  
  Iy = Iy + Iadjyr
  
  If (     TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iys,Ins,Ids,Ihs,Ims,Iss)) Go To 1
  
  If (.not.TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iye,Ine,Ide,Ihe,Ime,Ise)) Then
     Write (LM,'("End time found:",1x,i2.2,"/",i2.2,"/",i4,1x,3i2.2)') In, Id, Iy, Ih, Im, Is
     Return
  End If
  
  T_Init = Timz(Ihi,Imi,Isi)
  If (Idi .gt. Isday) T_Init = T_Init + 86400.0
  Time  = Timz(Ih,Im,Is)
  If (Ids .gt. Isday) Time = Time + 86400.0
  
  Olat1 = Olat
  Olon1 = Olon
  
  !  Move the grid if needed
  
  If (Nmosm .eq. 1) Then
     Olat1 = Olat + Sv*(Time-T_Init)/ 111120.0
     Olon1 = Olon + Su*(Time-T_Init)/(111120.0*Cos(Olat*ctr))
  End If
  
  !  Position of radar [km] relative to the origin
  
  Rx_Pos = (Rlon - Olon1) * 111.12 * Cos(Rlat*Ctr)
  Ry_Pos = (Rlat - Olat1) * 111.12
  
  !  Print out some info at the beginning of each sweep
  
  If (Isn .ne. Isn_old) Then
     Call Ctme(Time,Ih,Im,Is)
     hfrz = -999.0
     If (Ivpr .eq. 1) Call Get_Freezing_Level(Iy,In,Id,Ih,Im,Is,hfrz)
     
     Write (LM,10) Iy,In,Id, Ih, Im, Is, Azm, Elev, Rx_Pos, Ry_Pos,  &
          Olat1,Olon1, Isn, Irn, Alt, Ialt, Rdel, Gtln, Xnyq, N_Gates, Rname, Ntotal_Azms, hfrz
     
10   Format (/I4.4,"/",2i2.2,1x,3i2.2,' Az:',f5.1,' El:',f5.2, &
          1x,'Xz,Yz:',2f7.1,' Corner Pt:',2f9.3,' Sweep #:',i3, 1x,'Rec #:',i6, &
          1x,'Alt:',F7.1,' Ialt:',i6,' Rdel:',F7.3,' Gtln:',F7.3, &
          1x,'Nyq:',F7.2,' Gates:',i5,1x,a, 1x,'Total Azms in previous sweep: ',i4, &
          1x,'Freezing level:',F8.3)
     Isn_old = Isn
     Ntotal_Azms = 0
  End If

  Iazm = Nint(Azm*Azm_Step) + 1
  Iazm = Mod(Iazm,360) + 1
  
  If (Iazm .lt. 0 .or. Iazm .gt. 360) Then
     Write (LM,'("Azimuth strange: ",I5)') Iazm
     Call Exit
  End If
  
  !  Find any sectors that require special elevations
  
  Nstep_Offset = -999
  Nstep_Max = -999
  
  If (Ioffset .gt. 0) Then
     Do i = 1, Ioffset
        If (Azm .ge. Sectors(1,i) .and. Azm .le. Sectors(2,i)) Then
           Nstep_Offset = NOffsets(i)
           Nstep_Max = MaxS(i)
        End If
     End Do
  End If
  
  Ntotal_Azms = Ntotal_Azms + 1
  
  Do Igate = 1, N_Gates
     IHy_step = Hybrid(Iazm,Igate)
     
     If (Ihy_step .eq. -999) Cycle     ! Don't do the rest of the loop
     
     If (Nstep_Offset .ne. -999) Then
        IHy_step = Ihy_step + Nstep_offset
        If (IHy_step .gt. Nstep_Max) IHy_step = Nstep_Max
     End If
     
     !     Check to make sure the elevation angle desired is within range of
     !     the VCP that is being used
     
     If (IHy_step .gt. 14) IHy_step = 14
     If (IHy_step .lt. 1) IHy_step = 1
     
     Elev_Min = Float(nVCP12(IHy_step))/10.0
     dbz = Ref(Igate)
     
     If (Idbz .eq. 1) dbz = ZT(Igate)
     
     Eldif = Elev - Elev_Min
     
     If (Eldif .ge. -0.1 .and. Polar(Iazm,Igate,1) .eq. -998.0) Then 
        If (dbz .ne. -999.0) dbz = dbz + dbz_offset
        
        If (Ivpr .eq. 1) Then
           Range = Float(Igate)*Gtln + Rdel
           
           Call Get_Topo(Azm, Elev, Range, Rx_Pos, Ry_Pos, Sx, Sy, Topo, Imax, Jmax, Hgt_MSL)
           
           dbz_correction = 0.0
           dbz_raw = dbz
           
           If (dbz .ne. -999.0) Then
              Call Sergey_VPR(dbz, dbzc, Elev, Range, Alt, hfrz, Hgt_MSL,0)
              dbz_correction = dbzc-dbz
              total_vpr_corr = total_vpr_corr + dbz_correction
              ntotal_vpr = ntotal_vpr + 1
              dbz = dbzc
           End If
        End If
        
        If (dbz .ne. -999.0) Then
           Polar(Iazm,Igate,1) = dbz
           Polar(Iazm,Igate,2) = elev
           Polar(Iazm,Igate,3) = DR(Igate)
           Polar(Iazm,Igate,4) = KD(Igate)
        End If
     End If
  End Do
  
  Go To 1
  
End Subroutine Hreflec

Subroutine Sergey_VPR(dbz, dbzc, Elev, Range, Alt, hfrz, Hgt_MSL, Iprt)
  
  !--------------------------------------------------------------
  !     Routine to implement the Sergey Matrosov VPR correction.
  !     Input variables:
  !     dbz - original reflectivity observation [dBZ]
  !     dbzc - corrected reflectivity [dBZ]
  !     Elev - elevation angle [deg from horizontal]
  !     Range - Range from radar
  !     Alt - Altitude of radar
  !     hfrz = Height of the freezing level
  !     Hgt_MSL = Topo height
  !     Iprt is a flag to trigger some diagnostic print outs
  
  !     DelZ is the bright band enhancement [dBZ]
  !     bbthickness is the depth of the bright band [km]
  !--------------------------------------------------------------
  
  Real, Parameter :: dbz_max_correction= 12.0
  Real, Parameter :: DelZ = 9.0
  Real, Parameter :: BBthickness = 0.700
  
  Common /VRP/ DelZ1, tanb
  
  !  Is freezing level defined? If not, return
  
  If (hfrz .eq. -999.0) Then
     dbzc = dbz
     Return
  End If
  
  !     Is the topo height defined? If not the point is outside the grid so return
  
  If (Hgt_MSL .eq. -999.0) Then
     dbzc = dbz
     Return
  End If
  
  !      DelZ = 6.8 - 0.05*Range
  !      If (DelZ .lt. 0.0) DelZ = 0.0
  !     h0 is the heigth of the bottom of the bright band [km]
  !     htop is the height of the top of the bright band [km]
  !     hmid is the height of the middle of the bright band [km]
  !     Hgt is the beam height [km]
  !     Max_Correction is the maximum dbz correction
  
  htop = hfrz + 0.200 
  Hgt = Beam_Hgt(Range, Elev) + Alt/1000.0
  h0 = htop - bbthickness
  If (h0 .lt. 0.0) h0 = 0.0
  
  hmid = htop - (bbthickness/2.0)
  If (hmid .lt. 0.0) hmid = 0.0
  
  !  if the beam height is less than h0 then no correction
  
  If (hgt .lt. h0) Then
     dbzc = dbz
     Return
  End If
  
  !     If the beam height is between h0 and hmid then correct for bottom half of bright band 
  
  If (Hgt .gt. h0 .and. Hgt .le. hmid) Then
     dbzc = dbz - 2.0*DelZ * ((Hgt - h0)/(htop-h0))
     Return
  End If
  
  !     If the beam height is between hmid and htop then correct for top half of bright band
  
  If (Hgt .gt. hmid .and. Hgt .le. htop) Then
     dbzc = dbz - DelZ + 2.0*(DelZ+DelZ1) * ((Hgt-0.5*h0-0.5*htop)/(htop-h0))
     Return
  End If
  
  !     If the beam is above htop then correct for the "Ice Slope". If the topo is above the freezing level
  !     then correct only down to the height of the topo not htop.
  
  H1 = Hgt_MSL/1000.0
  
  If (H1 .le. htop) Then
     dbz_correction = Delz1 + tanb*(Hgt-htop)
  Else
     dbz_correction = tanb*(Hgt-H1)
  End If
  
  ! Keep the dbc correction under control if extrapolating over large height differences
  
  If (dbz_correction .gt. dbz_max_correction) dbz_correction = dbz_max_correction
  
  dbzc = dbz + dbz_correction
  
  If (Iprt .eq. 1) Then
     Write (6,*) H1, htop, Hgt, Hgt-htop, dbz, dbz_correction, dbzc       
  End If
  
  Return
End Subroutine Sergey_VPR

Subroutine Get_Topo(Azm, Elev, Range, Rx_Pos, Ry_Pos, Sx, Sy, Topo, Imax, Jmax, Hgt_MSL)
  
  !     Routine to find the MSL height of a given point in polar
  !     coordinates from the radar
  
  Data Ctr /0.0174532925/

  Real, Dimension(Imax,Jmax) :: Topo
  
  !   Compute the x coordinate (relative to the corner) [km]
  
  x = Rx_pos + Sin(Azm*Ctr) * Range
  
  !   Calculate the y coordinate (relative to the corner) [km]
 
  y = Ry_pos + Cos(Azm*Ctr) * Range
  
  i = Nint(x/Sx) + 1
  j = Nint(y/Sy) + 1
  
  If (i .lt. 1 .or. i .gt. Imax .or. j .lt. 1 .or. j .gt. Jmax) Then
     Hgt_MSL = -999.0
     i_alt = -999
     j_alt = -999
     Return
  End If
  
  Hgt_MSL = Topo(i,j)
  
  Return
End Subroutine Get_Topo

Subroutine RdTopo(TopoF, Topo)
  
  !  Routine to read in the topographic gridded data
  
  Common /One/ Imax, Jmax
  Common /Two/ Sx, Sy, Olat, Olon, Nmosm, Su, Sv, Xnyq, Alt, dbz_offset, Ivpr, Iatten
  Common /Write/ LM
  
  Dimension Topo(Imax,Jmax)
  Integer FirstChar
  
  Character*100 TopoF, Line, File, FileF
  
  nch = FirstChar(TopoF,' ') - 1
  File = TopoF(1:nch)
  
  Open (66, File=File,Iostat=Ierr,Err=1010,Status='Old')
  
  Read (66,'(a)') Line
  nch = FirstChar(Line,' ') - 1
  
  If (Line(1:9) .ne. 'topoe.prm') Then
     Write (LM,'(/"Invalid topo .prm file",1x,a)') Line(1:nch)
     Call Exit
  End If
  
  Read (66,'(a)') FileF
  
  ncht = FirstChar(FileF,' ') - 1
  File = FileF(1:ncht)
  
  Read (66,*) Tlat, Tlon
  Read (66,*) Tsx, Tsy
  Read (66,*) Nx, Ny
  
  Write (LM,'(/" Topography File Name: ",a)') FileF(1:ncht)
  
  Write (LM,'(/" Topo Grid Info:"/ &
       10x,"Lower Left Corner:",1x,2F9.3/ &
       10x,"Resolution (m):",1x,2F7.1/ &
       10x,"Grid Size (Pts,Km):",1x,2I5,1x,"(",2F6.1,")")') &
       Tlat, Tlon, Tsx*1000, Tsy*1000, Nx, Ny, Nx*Tsx, Ny*Tsy
  
  Close (66)
  
  If (Nx .ne. Imax .or. Ny .ne. Jmax) Then
     Write (LM,'(" Topo & Radar grids do not match",4I6)') Nx, Imax, Ny, Jmax
     Call Exit
  End If
  
  If (Tsx .ne. Sx .or. Tsy .ne. Sy) Then
     Write (LM,'(" Topo & Radar resolutions do not match",4F6.1)') Tsx, Sx, Tsy, Sy
     Call Exit
  End If
  
  If (Tlat .ne. Olat .or. Tlon .ne. Olon) Then
     Write (LM,'(" Topo & Radar grid corners do not match",4F9.3)') Tlat, Olat, Tlon, Olon
     Call Exit
  End If
  
  Lrec = Nx * 4
  
  !  topo files are always written "big_endian" because the program that 
  !  creates the topo.ele files (topoe.c) writes them that way
  
  Open (66,Err=1010, File=File, Iostat=Ierr, Status='Old', Access='Direct', Recl=Lrec, Convert='Big_Endian')
  Nrec = 1
  Tot = 0.0
  Knt = 0
  
  Do j = 1, Ny
     Read (66,Iostat=Ierr,Rec=Nrec) (Topo(i,j),i=1,Nx)
     Nrec = Nrec + 1
     
     Do i = 1, Nx
        Knt = Knt + 1
        Tot = Tot + Topo(i,j)
     End Do
  End Do
  
  AveTopo = Tot/Float(Knt)
  
  Write (Lm,'("Number of topo records read:",i6," Ave Topo hgt:",F8.1/)') Knt, AveTopo
  
  Close (66)
  Return
  
1010 Write (6,'("Error:",i5," on file:",a)') Ierr, File
  Call Exit
End Subroutine RdTopo

Subroutine Get_Freezing_Level (Iy,In,Id,Ih,Im,Is,hfrz)
  
  !--------------------------------------------------------------
  !     Routine to return the height of the freezing level
  !     Interpolation is performed if radar time is between two freezing
  !     level times
  !--------------------------------------------------------------
  
  Common /Freezing/ Fr_Hgts(1000,2), Nfrz, Freezing_File
  Common /Write/ LM
  
  Character*100 Freezing_File
  
  Logical*4 TimeBefore
  
  hfrz = -999.0
  i = 0
  
1 i = i + 1
  
  If (i .ge. Nfrz) Then
     hfrz = Fr_Hgts(Nfrz,2)
     Return
  End If
  
  !--------------------------------------------------------------
  !     Convert the radar time (yymmdd hhmmss) to fractional Julian Dates
  !     to match the time recorded with the freezing level heights
  !--------------------------------------------------------------
  
  JD = JulDy(In,Id,Iy)
  JD = JD + 365*(Iy-2005)
  Secs = Timz(Ih,Im,Is)
  Day = Float(JD) + Secs/86400.0
  
  DayJ2 = Fr_Hgts(i,1)
  
  !  2005 or 2006?
  
  Iyear = 2006
  If (DayJ2 .gt. 320.0) Iyear = 2005
  DayJ2 = DayJ2 + 365.0*(Float(Iyear-2005))
  If (DayJ2 .lt. Day) Go To 1
  
  If (i .eq. 1) Then
     hfrz = Fr_Hgts(1,2)
     Return
  End If
  
  hfrz1 = Fr_hgts(i-1,2)
  hfrz2 = Fr_hgts(i,2)
  
  DayJ1 = Fr_Hgts(i-1,1)
  Iyear = 2006
  If (DayJ1 .gt. 320.0) Iyear = 2005
  DayJ1 = DayJ1 + 365.0*(Float(Iyear-2005))
  
  !  Interpolate the freezing level to the radar time
  
  Fac = (Day-DayJ1)/(DayJ2-DayJ1)
  hfrz = hfrz1 + (hfrz2-hfrz1) * Fac
  
  Return
End Subroutine Get_Freezing_Level

Subroutine Load_Freezing_Level
  
  Common /Freezing/ Fr_Hgts(1000,2), Nfrz, Freezing_File
  Common /Write/ LM
  Common /Proj/ Project
  
  Character*16 Project
  Character*100 File, Line, Freezing_File
  
  Integer FirstChar
  
  !      File = '/Users/davej/projects/hmt/iop7/freezing-level-estimates-iop-07.asc'
  
  !      File =
  !     $   '/Users/davej/projects/hmt/hmt_data/freezing_heights_ruc.txt'
  
  File = Freezing_File
  
  !      If (Project(1:9) .eq. 'HMT-06/07') Then
  !         File =
  !     $   '/Users/davej/projects/hmt06/freezing_level_RUC_2006.txt'
  !      End If
  
  Open (45, File=File, Err=1010, Iostat=Ierr, Status='Old')
  
  Write (LM,'(/" Freezing level file: ",a)') File
  
  !  Skip the first header line
  
  Read (45,'(a)') Line
  
  i = 0
100 Read (45,'(a)',End=200) Line
  If (Line(1:2) .eq. '  ') Go To 100
  
  Read (Line,*) DayJ, F_lvl
  
  i = i + 1
  
  If (i .gt. 1000) Then
     Write (6,'("Too many freezing level hgts . aborting")')
     Call Exit
  End If
  
  Fr_Hgts(i,1) = DayJ
  Fr_Hgts(i,2) = F_lvl
  
  Go To 100
  
200 Nfrz = i
  
  Write (LM,'("Number of freezing level times:",i5/)') nfrz
  Return
  
1010 Write (6,'("Error:",i5," on file:",a)') Ierr, File
  Call Exit
  
End Subroutine Load_Freezing_Level

Subroutine Readuf_File (Azm, Elev, VR, DZ, ZT, DR, KD, Max_Gates, N_Gates, Iy, In, Id, Ih, Im, Is, Rdel, &
     Gtln, Rlat, Rlon, Ialt, Isn, Irn, Rname, EOF)
  
  !  Subroutine to read UF format tapes are return VE and DZ
  !  Variables of the subroutine:
   
  !      Man_Header  - Mandantory Header  Length = 45 words
  !      Iopt_Header - Optional Header    Length = Lenoh (max of 20)
  !      Loc_Header  - Local Header       Length = Lenlh (max of 20)
  !      Ifld_Header - Field Header       Length = Lenfh (max of 30)
  !      Idat_Header - Data Header        Length = Lendh (max of 17)
  
  !  Seven fields are currently supported for the SMART-R data . . . . . . .
  !         Name                Field
  !         ----                ------
  !          VE or VR     Doppler radial velocity [m/s]
  !          SW           Spectral width [m/s]
  !          DZ           Reflectivity [dBZ] with ground clutter cancellation
  !          ZT           Raw Reflectivity [dBZ]
  !          SQ           Signal Quality Index [0-1]
  !          DR           Differential Reflectivity [dBZ]
  !          KD           Kdp
  !          PH           PhiDP
  !          RH           RhoDP
  
  Integer, Parameter :: Lrec = 13500, MaxF = 7, Lenoh_max = 20, Lenlh_max = 20, Lendh_max = 17, Lenfh_max = 40
  
  Character(Len=8) Rname, SitNam, GenNam
  Character(Len=2) :: Char2, Ifname(MaxF)
  
  Logical :: EOF
  
  Real, Dimension(Max_Gates) :: DZ, VR, ZT, DR, KD, PH, dbz_corr, zdr_corr, CDZ, CDR
  
  Integer Man_Header(45), Iopt_Header(Lenoh_max), Loc_Header(Lenlh_max), Ifld_Header(Lenfh_max,MaxF)
  Integer Idat(Max_Gates,MaxF), Lenfh(MaxF), Idat_Header(Lendh_max)
  Integer*2 Iobuf(Lrec)

  Common /Two/ Sx, Sy, Olat, Olon, Nmosm, Su, Sv, Xnyq, Alt, dbz_offset, Ivpr, Iatten  
  Common /For/ Ioffset, Sectors(2,25), NOffsets(25), MaxS(25), Idbz
  Common /Write/ LM
  Common /CPU/ Machine
  
  ! ----------------------------------------------------------------
  ! Initalize data fields to missing in case the UF fields are missing
  ! ----------------------------------------------------------------
  
  DZ = -999.0
  VR = -999.0
  ZT = -999.0
  DR = -999.0
  KD = -999.0
  
  EOF = .False.
    
50 Read (50,End=2,Err=3000,Iostat=Ierr) Char2, Iobuf(2), (Iobuf(k),k=3,Iobuf(2)) 
  nrec = nrec + 1
  
  !  Is this a Universal Format Tape?
  
  If (Char2 .ne. 'UF') Then
     Write (6,'("Tape not Universal Format")')
     Call Exit
  End If
  
  !  .....CHECK THAT SPECIFIED RECORD LENGTH IS LESS THAN DIM(Iobuf)=Lrec
  
  If (Iobuf(2) .gt. Lrec) Then
     Write (6,12) Iobuf(2), Lrec
12   Format ('Record length=',I9,' exceeds buffer size=',I5)
     Call Exit
  End If
  
  !  .....Assign mandatory header (45 words)
  
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
  If (Azm .lt. 0.0) Azm = Azm + 360.0      ! Correct for the wierd way IRIS makes UF files
  
  Iy = Man_Header(26)                      ! Year
  In = Man_Header(27)                      ! Month
  Id = Man_Header(28)                      ! Day
  Ih = Man_Header(29)                      ! Hour
  Im = Man_Header(30)                      ! Minute
  Is = Man_Header(31)                      ! Second
  
  !  A good read? Check the date/time
  
  If (Iy .lt. 0 .or. In .lt. 0 .or. Id .lt. 0) Go To 50
  If (Ih .lt. 0 .or. Ih .gt. 23) Go To 50
  
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
  
  If (Machine .eq. 1) Then   ! switch the bits around if an intel Mac
     Do i = 11, 14
        Iobuf(i) = IShftC(Iobuf(i), 8, 16)
     End Do
  End If
  
  Write (Rname,'(4a2)') (Iobuf(i),i=11,14)
  Call fchar(Rname)           ! get rid of non-characters in the string

  !  Site Name . . . . . . . . . . . . .
  
  If (Machine .eq. 1) Then   ! switch the bits around if an intel Mac
     Do i = 15, 18
        Iobuf(i) = IShftC(Iobuf(i), 8, 16)
     End Do
  End If
  
  Write (SitNam,'(4a2)') (Iobuf(i),i=15,18)
  Call fchar(SitNam)           ! get rid of non-characters in the string

  !  Generating Facility Name . . . . . 
  
  If (Machine .eq. 1) Then    ! switch the bits around if an intel Mac
     Do i = 41, 44
        Iobuf(i) = IShftC(Iobuf(i), 8, 16)
     End Do
  End If
  
  Write (GenNam,'(4a2)') (Iobuf(i),i=41,44)
  Call fchar(GenNam)          ! get rid of non-characters in the string
  
  !  Correct for Y2K problems with NCAR translators
  
  Iy2k = 100
  
  If (GenNam(1:4) .eq. 'IRIS') Iy2k = 2000
  Iy  = Iy  - Iy2k
  Igy = Igy - Iy2k      
  
  !  Radar Latitude . . . . . . . . . . . . .
  
  Rlat = Float(Man_Header(19)) + Float(Man_Header(20))/60.0 + Float(Man_Header(21))/230400.0
  
  !  Radar Longitude . . . . . . . . . . . . . . .
  
  Rlon = Float(Man_Header(22)) + Float(Man_Header(23))/60.0 + Float(Man_Header(24))/230400.0
  Rlon = -Abs(Rlon)
  
  !  Radar altitude above sea level
  
  Ialt = Man_Header(25)
  
  !      Write (Char2,'(a2)') Man_Header(32)
  
  ! .....Assign Optional Header
  
  Lenoh = Man_Header(4) - Man_Header(3)
  
  If (Lenoh .gt. Lenoh_max) Then
     Write (6,'("Lengh of optional header too long:",i5)') Lenoh
     Call Exit
  End If
  
  If (Lenoh .gt. 0) Then
     Do i=1,Lenoh
        Iopt_Header(i) = Iobuf(45+i)
     End Do
  End If
  
  ! .....Assign Local Use header info
  
  Jstart = Man_Header(4) - 1
  Lenlh = Man_Header(5) - Man_Header(4)
  
  If (Lenlh .gt. Lenlh_max) Then
     Write (6,'("Lengh of local header too long:",i5)') lenlh
     Call Exit
  End If
  
  If (Lenlh .gt. 0) Then
     Do i=1, Lenlh
        Loc_Header(i) = Iobuf(Jstart+i)
     End Do
  End If
  
  ! .....Assign data header info
  
  Jstart = Man_Header(5)
  Nfe = Iobuf(Jstart)       ! Number of data fields
  
  If (Nfe .gt. MaxF) Then
     Write (6,'(i2," data fields currently allowed - need ",i5)') MaxF, Nfe 
     Call Exit
  End If
  
  Lendh = 3 + 2*Nfe
  
  If (Lendh .gt. Lendh_max) Then
     Write (6,'("Length of data header too long:",i5)') Lendh
     Call Exit
  End If
  
  Do j=1, Lendh
     Idat_Header(j) = Iobuf(Jstart + j - 1)
  End Do
  
  ! .....ASSIGN FIELD NAME AND FIELD HEADER INFO
  
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
     !        GET FIELD NAME
     
     If (Machine .eq. 1) Then  ! switch the bits around if an intel Mac
        Idat_Header(2+(2*J)) = IShftC(Idat_Header(2+(2*J)), 8, 16)
     End If
     
     Write (Ifname(j),'(a2)') Idat_Header(2+(2*J))
     
     !        GET STARTING WORD AND LENGTH FOR THIS FIELD
     
     Jstart = Idat_Header(3+(2*J))
     Lenfh(J) = Iobuf(Jstart) - Jstart
     Jstart = Jstart - 1
     
     If (Lenfh(j) .gt. Lenfh_max) Then
        Write (6,'("Length of field header too big:",i5)') Lenfh(j)
        Call Exit
     End If
     
     !   Define Field Header
     
     Do I = 1, Lenfh(j)
        Ifld_Header(i,j) = Iobuf(Jstart+I)
     End Do
  End Do
  
  !     GET STARTING POINTER FOR DATA IN EACH FIELD
  
  Do j = 1, Nfe
     Jstart = Ifld_Header(1,j) - 1
     N_Gates = Ifld_Header(6,j)
     
     !  Load up only as many gates as the main program has allocated for
     
     If (N_Gates .gt. Max_Gates) Then
        N_Gates = Max_Gates
     End If
     
     Do i = 1, N_Gates
        Idat(i,j) = Iobuf(Jstart+i)
     End Do
  End Do

  AdjZ = -999.0
  Gtz = -999.0
  N_GatesZ = -999
  RdelZ = -999.0
  
  Do j = 1, Nfe     
     If (Ifname(j) .eq. 'ZT' .or. Ifname(j) .eq. 'RE') Then    ! total Reflectivity
        Scale = Ifld_Header(2,j)
        RdelZ = Float(Ifld_Header(3,j))
        GtZ   = Float(Ifld_Header(5,j)) / 1000.0
        AdjZ  = Float(Ifld_Header(4,j)) / 1000.0
        N_GatesZ = Ifld_Header(6,j)
        RangeOne = RdelZ + AdjZ        ! range of gate one
        ij = 0
        If (N_GatesZ .gt. Max_Gates) N_GatesZ = Max_Gates  
        IF (N_GatesZ .le. 0) N_GatesZ = 1
        
        Do i = 1, N_GatesZ
           Range= Float(i) * GtZ + RangeOne
           
           If (Range .gt. 0.0) Then     ! only use positive ranges
              ij = ij + 1
              ZT(ij) = Float(Idat(i,j))/Scale
              If (Idat(i,j) .eq. Ibad) ZT(ij) = -999.0
           End If
        End Do
        
        !  Readjust some values in case first gate(s) discarded because of (-) range
        
        RdelZ = RdelZ + Float(N_GatesZ - ij) * GtZ ! Readjust RdelR
        N_GatesZ = ij                              ! Readjust N_gatesR
     End If   ! For ZT
     
     If (Ifname(j) .eq. 'DZ' .or. Ifname(j) .eq. 'RE') Then     !  the filtered Reflectivity field
        Scale = Ifld_Header(2,j)
        RdelR = Float(Ifld_Header(3,j))
        GtR   = Float(Ifld_Header(5,j))/1000.0
        AdjR  = Float(Ifld_Header(4,j))/1000.0
        N_GatesR = Ifld_Header(6,j)
        RangeOne = RdelR + AdjR        ! range of gate one
        ij = 0
        If (N_GatesR .gt. Max_Gates) N_GatesR = Max_Gates
        IF (N_GatesR .le. 0) N_GatesR = 1
        
        Do i = 1, N_GatesR
           Range = Float(i) * GtR + RangeOne
           
           If (Range .ge. 0.0) Then     ! only use positive ranges
              ij = ij + 1
              DZ(ij) = Float(Idat(i,j))/Scale
              If (Idat(i,j) .eq. Ibad) DZ(ij) = -999.0
           End If
        End Do
        
        !  Readjust some values in case first gate(s) discarded because of (-) range
        
        RdelR = RdelR + Float(N_GatesR - ij) * GtR ! Readjust RdelR
        N_GatesR = ij       ! Readjust N_GatesR
        
     End If   ! for DZ
     
     If (Ifname(j) .eq. 'DR') Then     !  the 'DR' field
        Scale = Ifld_Header(2,j)
        RdelR = Float(Ifld_Header(3,j))
        GtR   = Float(Ifld_Header(5,j))/1000.0
        AdjR  = Float(Ifld_Header(4,j))/1000.0
        N_GatesR = Ifld_Header(6,j)
        RangeOne = RdelR + AdjR        ! range of gate one
        ij = 0
        If (N_GatesR .gt. Max_Gates) N_GatesR = Max_Gates
        IF (N_GatesR .le. 0) N_GatesR = 1
        
        Do i = 1, N_GatesR
           Range = Float(i) * GtR + RangeOne
           
           If (Range .ge. 0.0) Then     ! only use positive ranges
              ij = ij + 1
              DR(ij) = Float(Idat(i,j))/Scale
              If (Idat(i,j) .eq. Ibad) DR(ij) = -999.0
           End If
        End Do
        
        !  Readjust some values in case first gate(s) discarded because of (-) range
        
        RdelR = RdelR + Float(N_GatesR - ij) * GtR ! Readjust RdelR
        N_GatesR = ij       ! Readjust N_GatesR
        
     End If   ! for DR
     
     If (Ifname(j) .eq. 'KD') Then     !  the 'KD' field
        Scale = Ifld_Header(2,j)
        RdelK = Float(Ifld_Header(3,j))
        GtK   = Float(Ifld_Header(5,j))/1000.0
        AdjK  = Float(Ifld_Header(4,j))/1000.0
        N_GatesK = Ifld_Header(6,j)
        RangeOne = RdelK + AdjK        ! range of gate one
        ij = 0
        If (N_GatesK .gt. Max_Gates) N_GatesK = Max_Gates
        IF (N_GatesK .le. 0) N_GatesK = 1
        
        Do i = 1, N_GatesK
           Range = Float(i) * GtK + RangeOne
           
           If (Range .ge. 0.0) Then     ! only use positive ranges
              ij = ij + 1
              KD(ij) = Float(Idat(i,j))/Scale
              If (Idat(i,j) .eq. Ibad) KD(ij) = -999.0
           End If
        End Do
        
        !  Readjust some values in case first gate(s) discarded because of (-) range
        
        RdelK = RdelK + Float(N_GatesK - ij) * GtK ! Readjust RdelR
        N_GatesK = ij       ! Readjust N_GatesK
        
     End If   ! for KD
     
     !  Get Velocity data
     
     If (Ifname(j) .eq. 'VR' .or. Ifname(j) .eq. 'VE') Then      ! the 'VR' field
        Scale = Ifld_Header(2,j)
        RdelV = Float(Ifld_Header(3,j))
        IgtV   = Ifld_Header(5,j)
        GtV = Float(IgtV)/1000.0
        AdjV  = Float(Ifld_Header(4,j))/1000.0
        N_gatesV = Ifld_Header(6,j)
        N_GatesVOrig = N_GatesV
        Xnyq = Float(Ifld_Header(20,j))/Scale
        RangeOne = RdelV + AdjV        ! range of gate one
        ij = 0
        If (N_GatesV .gt. Max_Gates) N_GatesV = Max_Gates
        IF (N_GatesV .le. 0) N_GatesV = 1
        
        Do i = 1, N_gatesV
           Range = Float(i)* GtV + RangeOne
           
           If (Range .ge. 0.0) Then
              ij = ij+1
              VR(ij) = Float(Idat(i,j))/Scale
              If (Idat(i,j) .eq. Ibad) VR(ij) = -999.99
           End If
        End Do               ! i=1,ngates
        
        !  Readjust some values in case first gate(s) discarded because of (-) range
        
        RdelV = RdelV + Float(N_GatesV -ij) * GtV ! Readjust RdelV
        N_GatesV = ij                             ! Readjust N_gatesV
     End If   ! for VR or VE

     If (Ifname(j) .eq. 'PH') Then   ! for PhiDP [deg]
        Scale = Ifld_Header(2,j)
        RdelP = Ifld_Header(3,j)
        GtP   = Float(Ifld_Header(5,j))/1000.0
        N_GatesP = Ifld_Header(6,j)
        
        Do i = 1, N_gatesP
           Spec = Idat(i,j)
           PH(i) = Abs(Spec/Scale)
           If (Idat(i,j) .eq. Ibad) PH(i) = -999.0
        End Do
     End If   ! for PH
  End Do                     ! j=1,Nfe

  If (Iatten .eq. 1) Then
     Call atten_corr (DZ, DR, PH, CDZ, CDR, dbz_corr, zdr_corr, Nsteps, N_GatesZ, -999.0)

     Do i = 1, N_GatesZ
        DZ(i) = CDZ(i)
        DR(i) = CDR(i)
     End Do
  End If
 
  Rdel = RdelZ + AdjZ
  Gtln = GtZ
  N_Gates = N_GatesZ
  
  Return
  
2 EOF = .True.
  
  Return
  
3000 Write (6,'("Error on UF read: ",i5)') Ierr
  
  Call Exit
  
End Subroutine Readuf_File

Subroutine Cnvt_Cartesian (Polar, Nazms, Nbins, Rdel, Gtln, Azm_Step, Cart, &
     Imax, Jmax, Sx, Sy, Olat, Olon, Rlat, Rlon, Hinf) 
  
  !--------------------------------------------------------------
  !     Routine to interpolate a polar coordinate grid to a Cartesian
  !     coordinate grid using Cressman interpolation with a variable
  !     influence radius with range.
  
  !     The horizontal influence angle is given by Hinf and as a function
  !     of range the horizontal radius of influence is 
  !     R_infh = Tan(Hinf*Ctr) * Range
  
  !     Input variables:
  !        Polar - original polar coordinate data (azms,bins,2) 1:dbz 2:elev
  !	 Nazms - number of azimuths
  !	 Nbins - number of range gates
  !	 Rdel - Range delay [km]
  !	 Gtln - gate spacing [km]
  !	 Azm_Step - azimuthal spacing between rays [deg]
  !	 Cart - Cartesian array
  !	 Imax - number of east-west (x axis) points
  !	 Jmax - number of north-south (y axis) points
  !	 Sx - x axis spacing [km]
  !	 Sy - y axis spacing [km]
  !	 Olat - Cartesian grid latitude origin (lower left corner) [deg]
  !	 Olon - Cartesian grid longitude origin (lower left corner) [deg]
  !	 Rlat - Radar latitude [deg]
  !	 Rlon - Radar longitude [deg]
  !	 Hing - horizontal influence radius (deg)
  !--------------------------------------------------------------

  ! There are 4 fields in Polar 1=dbz 2=elev 3=Zdr 4=Kdp
  ! There are 3 field in the Cartesian array 1=dbz 2=Zdr 3=Kdp

  Real, Dimension(Nazms,Nbins,4) :: Polar
  Real, Dimension(Imax,Jmax,3) :: Weight, Cart
  
  Data Ctr /0.0174532925/
  
  !-----------------------------------------------------------  
  ! Initialize the Cartesian & weighting fuction arrays
  !-----------------------------------------------------------

  Weight = 0.0
  Cart = -999.0

  !-----------------------------------------------------------   
  !  Position of radar [km] relative to the origin
  !-----------------------------------------------------------  
  
  Rx_Pos = (Rlon - Olon) * 111.12 * Cos(Rlat*Ctr)
  Ry_Pos = (Rlat - Olat) * 111.12
  
  !-----------------------------------------------------------   
  !  The big loop
  !-----------------------------------------------------------   
  
  Do Iazm = 1, Nazms
     Azm = Float(Iazm-1) * Azm_Step
     
     Do Igate = 1, Nbins
        If (Polar(Iazm, Igate, 1) .eq. -999.0 .or.  &
             Polar(Iazm, Igate, 2) .eq. -999.0) Cycle ! Missing dbz or elev
        
        Elev = Polar(Iazm, Igate, 2)
        Range = Rdel + Float(Igate) * Gtln
        
        If (Range .le. 0.0) Cycle      ! in case the Rdel is negative
        
        R_infh = Tan(Hinf*Ctr) * Range
        
        !     Keep the influence radius from getting ridiculously small close to
        !     the radar
        
        If (R_infh .lt. 0.2) R_infh = 0.200
        
        !   Compute the x coordinate (relative to the corner) [km]
        
        x = Rx_pos + Sin(Azm*Ctr) * Cos(Elev*Ctr) * Range
        
        !   Calculate the y coordinate (relative to the corner) [km]
        
        y = Ry_pos + Cos(Azm*Ctr) * Cos(Elev*Ctr) * Range
        
        !   To cut down on the amount of arithmetic, find the i grid points
        !   within the horizontal influence radius 
        
        X_min = x - R_infh/2.0
        X_max = x + R_infh/2.0
        I_min = Nint(X_min/Sx) + 1
        I_max = Nint(X_max/Sx) + 1
        
        !   Make sure the limits are within the grid
        
        If (I_min .lt. 1) I_min = 1
        If (I_max .gt. Imax) I_max = Imax
        If (I_min .gt. Imax .or. I_max .lt. 1) Cycle
        
        !   To cut down on the amount of arithmetic, find the j grid points
        !   within the horizontal influence radius 
        
        Y_min = y - R_infh/2.0
        Y_max = y + R_infh/2.0
        J_min = Nint(Y_min/Sy) + 1
        J_max = Nint(Y_max/Sy) + 1
        
        !   Make sure the limits are within the grid
        
        If (J_min .lt. 1) J_min = 1
        If (J_max .gt. Jmax) J_max = Jmax
        If (J_min .gt. Jmax .or. J_max .lt. 1) Cycle
        
        !   Determine weighting function for all points within the influence radius
        
        Do i = I_min, I_max
           X_pos = Float(i-1) * Sx
           Del_x = x - X_pos
           
           Do j = J_min, J_max
              Y_pos = Float(j-1) * Sy
              Del_y = y - Y_pos
              D_sqr = Del_x*Del_x + Del_y*Del_y
              R_sqr = R_infh * R_infh
              
              !   If the polar points are within the influence radius then use them
              
              If (R_sqr .gt. D_sqr) Then
                 W = (R_sqr - D_sqr) / (R_sqr + D_sqr)
                 Z = Cart(i,j,1)
                 Zdr = Cart(i,j,2)
                 Xdp = Cart(i,j,3)
                 Znew = Polar(Iazm,Igate,1)
                 Zdrnew = Polar(Iazm,Igate,3)
                 Xdpnew = Polar(Iazm,Igate,4)
                 
                 If (Znew .ne. -999.0) Then                  
                    Z = Z * Weight(i,j,1) + Znew * W 
                    Z = Z / (Weight(i,j,1) + W)
                    Cart(i,j,1) = Z
                    Weight(i,j,1) = Weight(i,j,1) + W 
                 End If

                 If (Zdrnew .ne. -999.0) Then                  
                    Zdr = Zdr * Weight(i,j,2) + Zdrnew * W 
                    Zdr = Zdr / (Weight(i,j,2) + W)
                    Cart(i,j,2) = Zdr
                    Weight(i,j,2) = Weight(i,j,2) + W 
                 End If

                 If (Xdpnew .ne. -999.0) Then                  
                    Xdp = Xdp * Weight(i,j,3) + Xdpnew * W 
                    Xdp = Xdp / (Weight(i,j,3) + W)
                    Cart(i,j,3) = Xdp
                    Weight(i,j,3) = Weight(i,j,3) + W 
                 End If
              End If
           End Do
        End Do
     End Do
  End Do
  
  Return
End Subroutine Cnvt_Cartesian

Block Data
   
   Common /One/ Imax, Jmax
   Common /Two/ Sx, Sy, Olat, Olon, Nmosm, Su, Sv, Xnyq, Alt, dbz_offset, Ivpr, Iatten
   Common /Thr/ Nazms, Nbins, Azm_Step, RdelH, GtlnH, RlatH, RlonH
   Common /For/ Ioffset, Sectors(2,25), NOffsets(25), MaxS(25), Idbz
   Common /Five/ Itime_Limits(12), Init_Time(6), Isday, Iadjyr
   Common /Write/ LM
   Common /Interp/ Hinf
   Common /VRP/ DelZ1, tanb
   Common /Proj/ Project
   
   Data LM /64/
   
   Character Project*16
End Block Data
