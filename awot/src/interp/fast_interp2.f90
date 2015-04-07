Program Fast_Interp2

  !  Program to interpolate fore/aft scanning date to 3-D grids  
  !     such as is produced by the SOLO translator
  !     The program must be run for each .for & .aft files

  !  All output is directed to a file @.log where @ is the same as the .prm

  !  The parameters used in the program:

  !  Namof- name of the disc file that will contain output data
  !  Namdf- name of the disc file that contains the unfolded MERGED data
  !  Imax - maximum points in X direction (East-West)
  !  Jmax - maximum points in Y direction (North-South)
  !  Kmax - maximum points is Z direction (vertical)
  !  Sx   - distance between points (km) in X direction
  !  Sy   - distance between points (km) in Y direction
  !  Sz   - distance between points (km) in Z direction
  !  Z0   - height above surface of 1st plane (km)
  !         (Note: Definition of "surface" is defined at load time as that
  !         of the ground or mean sea level, concomitant with use of radar
  !         altitude or pressure altitude, respectively.)
  !  Nmosm- flag for type of grid reference desired
  !       - =0 with respect to the ground; in this case,
  !          OLAT - latitude of lower left hand corner
  !          OLON - longitude of lower left hand corner
  !       - =1 with respect to the moving grid center; in that case,
  !          Olat = Olat + Sv*(Time-T_init)/111190.0
  !          OLon = Olon + Su*(Time-T_init)/(111190.0*Cos(Olat))
  !  Sazm - starting azimuth
  !  Eazm - ending azimuth
  !  IHS, IMS, ISS - starting time of 1st leg
  !  IHE, IME, ISE - ending time of 1st leg

  !  Summary of disk lu's used:
  !       10 - parameter file
  !       LM - log file
  !       66 - terrain file
  !       95 - OUTPUT file (.for or .aft)
  !       75 - MERGED data file

  Integer*2 Ihead(45)
  Integer*4 Iheadr(40)
  Integer*4 FirstChar

  Logical*4 OKTime

  Common /One/ Imax, Jmax, Kmax
  Common /Two/ Sx, Sy, Sz, Z0, Olat, Olon, Sazm(2), Eazm(2), Nmosm, Su, Sv, Xniq
  Common /Tre/ Flid, Project
  Common /Five/ Itime_Limits(12), Init_Time(6), Isday, Iadjyr
  Common /Six/ Ctr, Eloff, CantDrft, Rotcr, PitchCor, GsCor, RangeMax 
  Common /Topo/ Tlat, Tlon, Tsx, Tsy, Max_topox, Max_topoy, Itopo
  Common /Write/ LM

  Character*100 TopoFile, Line, Ifile, PrmFile, Namdf, Namof
  Character*24 Time_Current
  Character*3 For, Aft, Itype
  Character*8 Flid
  Character*16 Project

  Data For /'for'/, Aft /'aft'/

  !    Open the data file containing the parameter information.

  n = IargC()

  If (n .lt. 1) Then
     Write (6,'("Enter prm file name:",$)')
     Read (5,'(a)') PrmFile
  Else
     Call GetArg(1,PrmFile)
  End If

  !  Open the parameter file

  nch = LenTrim(PrmFile)
  Ifile = Prmfile(1:nch)

  Open (10, File=Ifile, Err=1010, Iostat=Ierr, Status='Old')

  Write (6,'(/" P-3 Doppler Interpolation program -  .prm File:",a)') Ifile(1:nch)

  !  Create and open the .log file with the same name as the .prm file

  nc = LastChar(Prmfile,'.')
  Prmfile = Prmfile(1:nc) // 'log'
  nch = LenTrim(PrmFile)
  Ifile = Prmfile(1:nch)

  If (LM .ne. 6) Then
     Open (LM, File=Ifile(1:nch),Iostat=Ierr, Err=1010)
  End If

  Write (LM,'("  .log File:",a)') Ifile(1:nch)

  Call FDATE(Time_Current)

  Write (LM,'("       Starting fast_interp2 ",a)') Time_Current

  !  Read parameters from the disc file
  !     The first line is the indicator for the type of .prm file (fast_interp2.prm)

  Read (10,'(a)') Line
  nch = FirstChar(Line,' ') - 1

  If (Line(1:nch) .ne. 'fast_interp2.prm') Then
     Write (LM,'(/"  Improper type of .prm file - must start with fast_interp2.prm ",a/)') Line(1:nch)
     Stop 'fast_interp2'
  End If

  !  Next, read the line containing the type (for or aft) and topo info

  Read (10,'(a3,1x,i1,1x,a)') Itype, Itopo, TopoFile

  If (Itype .ne. For .and. Itype .ne. Aft) Then
     Write (LM,'(/"File type=",a," not ",a," or ",a," ABORTING!!!")') Itype, For, Aft
     Call Exit
  End If

  !  Is there a topo file to help remove the ground clutter??

  If (Itopo .eq. 1) Then
     nch = FirstChar(TopoFile,' ') - 1
     Ifile = TopoFile(1:nch)
     Open (66, File=Ifile(1:nch), Iostat=Ierr, Err=1010)

     Write (LM,'(/" Topography .prm File Name: ",a)') Ifile(1:nch)

     Read (66,'(a)') Line
     nch = FirstChar(Line,' ') - 1

     If (Line(1:nch) .ne. 'topoe.prm') Then
        Write (LM,'(/"Invalid topo .prm file",1x,a)') Line(1:nch)
        Call Exit
     End If

     Read (66,'(a)') Ifile

     nch = FirstChar(Ifile,' ') - 1
     Ifile = Ifile(1:nch)

     Read (66,*) Tlat, Tlon
     Read (66,*) Tsx, Tsy
     Read (66,*) Max_topox, Max_topoy

     Close (66)

     Write (LM,'(/" Topography File Name: ",a)') Ifile(1:nch)

     Write (LM,'(/"Topo Grid Info:"/ &
             10x,"Lower Left Corner:",1x,2F9.3/ &
             10x,"Resolution (m):",1x,2F6.1/ &
             10x,"Grid Size (Pts,Km):",1x,2I5,1x,"(",2F6.1,")")') &
             Tlat, Tlon, Tsx*1000, Tsy*1000, Max_topox, Max_topoy, Max_topox*Tsx, Max_topoy*Tsy

     Lrec = Max_topox * 4

     Open (66,Err=1010, File=Ifile, Iostat=Ierr, Status='Old', Access='Direct', Recl=Lrec)
  End If

  !  Next, the file name of the input (cleaned up) data

  Read (10,'(a)') Line
  nch = FirstChar(Line,' ') - 1
  Namdf = Line(1:nch)
  nch = LenTrim(Namdf)

  Write (LM,'(/" .cor file name: ",a)') Namdf(1:nch)

  !  Next, the generic name of the output files
  !  They will be named .for and .aft for the forward and aft scans,
  !  respectively

  Read (10,'(a)') Line
  nch = FirstChar(Line,' ') - 1
  Namof = Line(1:nch)
  nch = LenTrim(Namof)

  Write (LM,'(/" .for & .aft file names: ",a)') Namof(1:nch)

  !  The grid parameters

  Read (10,*) Imax, Jmax, Kmax, Sx, Sy, Sz, Z0
  Read (10,'(i1,1x,f6.0,1x,f7.0,1x,3i2,1x,3i2,1x,f5.0,1x,f5.0)') Nmosm, Olat, Olon, Iyi, Ini, Idi, Ihi, Imi, Isi, Su, Sv

  !  The azimuthal limits

  Read (10,*) Sazm(1), Eazm(1), Sazm(2), Eazm(2)

  !  The beam pointing & platform attitude/motion corrections

  Read (10,*) Rotcr, Eloff, CantDrft, PitchCor, GsCor, RangeMax

  If (RangeMax .le. 0.0) RangeMax = 100.0

  !  Lastly, The starting and ending times

  Read (10,'(3i2,1x,3i2,1x,3i2,1x,3i2)') Iys, Ins, Ids, Ihs, Ims, Iss, Iye, Ine, Ide, Ihe, Ime, Ise
  Close(10)

  Iadjyr = 1900
  If (Iys .lt. 75) Iadjyr = 2000

  If (.not. OKTime(Iyi+Iadjyr,Ini,Idi,Ihi,Imi,Isi)) Then
     Write (LM,'("Time: ",i2.2,"/"i2.2,"/",i2.2,1x,3i2.2," is not OK")') Iyi, Ini, Idi, Ihi, Imi, Isi
     Stop
  End If

  Init_Time(1) = Iyi
  Init_Time(2) = Ini
  Init_Time(3) = Idi
  Init_Time(4) = Ihi
  Init_Time(5) = Imi
  Init_Time(6) = Isi

  !  Are the start times OK?

  If (.not. OKTime(Iys+Iadjyr,Ins,Ids,Ihs,Ims,Iss)) Then
     Write (LM,'("Time: ",i2.2,"/"i2.2,"/",i2.2,1x,3i2.2," is not OK")') Iys, Ins, Ids, Ihs, Ims, Iss
     Stop
  End If

  !  Are the end times OK?

  If (.not. OKTime(Iye+Iadjyr,Ine,Ide,Ihe,Ime,Ise)) Then
     Write (LM,'("Time: ",i2.2,"/"i2.2,"/",i2.2,1x,3i2.2," is not OK")') Iye, Ine, Ide, Ihe, Ime, Ise
     Stop
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

  !     Open the .cor data file

  nch = LenTrim(Namdf)
  Ifile = Namdf(1:nch)

  Open (75,File=Ifile(1:nch), Form='UNFORMATTED', Err=1020, Iostat=Ierr, Status='Old')

  !     Read the first record to get the start time of the file
  !     and other parameters needed for the header file

1 Read (75,End=1020,Err=1020,Iostat=Ierr) (Ihead(i),i=1,45),(Iheadr(i),i=1,40)

  IsDay = Iheadr(3)

  !     Project and flight identifiers

  Write (Flid,'(4a2)') (Ihead(kk),kk=24,27)
  Write (Project,'(6a2,4H- - )') (Ihead(kk),kk=28,33)

  !     Nyquist interval (m/s)

  Xniq = Float(Ihead(17)) / 10.0
  Alt_Flag = Ihead (44)

  Rewind 75

  Xsiz = Sx * Float(Imax)
  Ysiz = Sy * Float(Jmax)

  Write (LM,32) Imax, Jmax, Kmax, Sx, Sy, Sz, Z0, Nmosm, Olat, Olon, Iyi, Ini, Idi, Ihi, Imi, Isi, &
       Su, Sv, Sazm(1), Eazm(1), Sazm(2), Eazm(2), Ihs, Ims, Iss, Ihe, Ime, Ise, Namof, Namdf, Xsiz, &
       Ysiz, Rotcr, Eloff, CantDrft, PitchCor, GsCor, Flid, Project, Xniq, RangeMax, Alt_Flag

32 Format (/'.cor file headers:'/ &
        ' Imax, Jmax, Kmax:',3I4/ &
        ' Sx Sy Sz Z0:',4F7.2/    &
        ' Nmosm Olat Olon:',I1,1X,2F8.3/ &
        ' Init Time: ', 3i2.2,1x,3i2.2/  &
        ' Storm u,v components: ',2f5.1/ &
        ' Start Azm, End Azm:',4F6.1/    &
        ' Start time - End time:',1X,3I2.2,1X,3I2.2/ &
        ' Output file dir:',1X,a/ &
        ' Input file name:',1X,a/ &
        ' Domain size (km):',1x,2f7.1/ &
        ' Rotation correction [deg]:',1x,f7.2/ &
        ' Tilt angle correction [deg]:',1x,f7.2/ &
        ' Canting+Drift angle offset [deg]:',1x,f7.2/ &
        ' Pitch angle offset [deg]:',1x,f7.2/ &
        ' Ground speed offset [m/s]:',1x,f7.2/ &
        ' Flight ID:',a/ &
        ' Project:',a/ &
        ' Nyquist Interval [m/s]:',F7.2/ &
        ' Max Range:',F7.1/ &
        ' Altitude Flag:',F5.0/)

  ! Write data to disc file

  nch = LenTrim(Namof)
  Namof = Namof(1:nch) // '.' // Itype
  nch = LenTrim(Namof)

  Ifile = Namof(1:nch)

  Open (95, Err=1020, File=Ifile(1:nch), Iostat=Ierr, Form='UNFORMATTED')

  Write (LM,'(/" Output File Name:",a/)') Namof(1:nch)

  !   Write the header data

  Rlat = -999.0
  Rlon = -999.0
  Dum = 0.0

  Call FDATE(Time_Current)

  Write (95) Time_Current, Imax, Jmax, Kmax, Sx, Sy, Sz, Flid, Olat, Olon, Z0, Itime_Limits, Nmosm, Init_Time, &
       Su, Sv, Project, Rotcr, Eloff, CantDrft, Namdf, Xniq, Rlat, Rlon, Itype, RangeMax, Alt_Flag, Dum, Dum

  Write (LM,12) Itype, Time_Current, Imax, Jmax, Kmax, Sx, Sy, Sz, Flid, Olat, Olon, Z0, Itime_Limits, Nmosm, Init_Time, &
       Su, Sv, Project, Rotcr, Eloff, CantDrft, Namdf, Xniq, RangeMax, Rlat, Rlon, Itype

12 Format (/'Header for .',a,' file:'/ &
        ' File created on: ',a/ &
        ' Imax Jmax Kmax: ',3i4/ &
        ' Sx Sy Sz [km]: ',3F7.3/ &
        ' Flight Id: ',a/ &
        ' Olat Olon Z0 :',2F8.3,F5.1/ &
        ' Start, End Times: ',3i2.2,1x,3i2.2,1x,3i2.2,1x,3i2.2/ &
        ' Nmosm:',I1,1x,'Init Time: ',3i2.2,1x,3i2.2/ &
        ' Storm u,v components: ', 2f7.2/ &
        ' Project: ',a/ &
        ' Azm Cor  El offset  CantDrft:',3f7.2/ &
        ' Name of .cor file: ',a/ &
        ' Nyquist: ',f7.2/ &
        ' Maximum Range Used: ',f7.2/ &
        ' Radar Position: ',2f8.2,' Type: ',a)

  Call Main_Body(Itype)

  Close (95)

  Call Exit

1010 Write (6,'("Error:",i5," on file:",a)') Ierr, Ifile
  Call Exit

1020 Write (LM,'("Error:",i5," on file:",a)') Ierr, Ifile
  Call Exit

End Program Fast_Interp2

Subroutine Main_Body(Itype)

  Common /One/ Imax, Jmax, Kmax
  Common /Write/ LM

  Dimension Data(7, Imax, Jmax, Kmax)

  Integer*2 Ihd1(Imax), Ihd2(Imax), Ihd3(Imax), Ihd4(Imax), Ihd5(Imax), Ihd6(Imax)

  Character*24 Time_Current
  Character*3 Itype

  !  Initialize the arrays that will hold the cartesian data

  Do l = 1, 7
     Do i = 1, Imax
        Do j = 1, Jmax
           Do k = 1, Kmax
              Data(l,i,j,k) = 0.0
           End Do
        End Do
     End Do
  End Do

  !  Do the 3-D interpolation

  Call FDATE(Time_Current)
  Write (LM,'(/" Starting interpolation ",a)') Time_Current

  Call Cint (Itype, Data)

  !  Close the input file

  Close (75)

  Call FDATE(Time_Current)
  Write (LM,'(/" Stopping interpolation ",a)') Time_Current
  Write (LM,'(/)')

  !  Write the interpolated data out to the disk

  Do k = 1, Kmax
     Write (LM,'(" Writing ",a," file: Plane #",i2)') Itype, k

     Do j = 1, Jmax
        Do i = 1, Imax
           Call GetData (i,j,k, Data, Trk, Vel, Azm, Ref, Rng, Tdif)
           Ihd1(i) = Nint(Trk * 10.0)  ! =ground-based azimuth
           Ihd2(i) = Nint(Azm * 10.0)  ! =tail azm "adjusted" for tilt
           Ihd3(i) = Nint(Vel * 10.0)  ! =radial velocity m/s
           Ihd4(i) = Nint(Ref * 10.0)  ! =reflectivity
           Ihd5(i) = Nint(Rng * 10.0)  ! =range to bin (km)
           Ihd6(i) = Nint(Tdif)        ! =time difference from init time
        End Do

        Write (95) (Ihd1(i), Ihd2(i), Ihd3(i), Ihd4(i), Ihd5(i), Ihd6(i),i=1,Imax)
     End Do
  End Do

  Return
End Subroutine Main_Body

Subroutine Cint (Itype, Data)

  Common /One/ Imax, Jmax, Kmax
  Common /Two/ Sx, Sy, Sz, Z0, Olat, Olon, Sazm(2), Eazm(2), Nmosm, Su, Sv, Xniq
  Common /Five/ Itime_Limits(12), Init_Time(6), Isday, Iadjyr
  Common /Six/ Ctr, Eloff, CantDrft, Rotcr, PitchCor, GsCor, RangeMax 
  Common /Topo/ Tlat, Tlon, Tsx, Tsy, Max_topox, Max_topoy, Itopo
  Common /Write/ LM

  Dimension Vel1(1024), Vel2(1024), Hed1(13), Hed2(13), Ref1(1024), Ref2(1024)
  Dimension Itime1(6), Itime2(6)

  Dimension Data(7,Imax,Jmax,Kmax)    
  Dimension Elev(Max_topox, Max_topoy)

  Character*3 Itype
  Logical*4 TimeOnBefore

  Equivalence (Ang1,Hed1(1)), (Plat1,Hed1(2)), (Plon1,Hed1(3)), (Gtln1,Hed1(4)), (Zgat1,Hed1(5)), (Azraw1,Hed1(6)),  &
       (SweepNum1,Hed1(7)), (Azm1,Hed1(8)), (Tilt1,Hed1(9)), (Alt1,Hed1(10)), (Rdel1,Hed1(11)), (Alt_Sel1,Hed1(12)), &
       (Ra1,Hed1(13))

  Equivalence (Ang2,Hed2(1)), (Plat2,Hed2(2)), (Plon2,Hed2(3)), (Gtln2,Hed2(4)), (Zgat2,Hed2(5)), (Azraw2,Hed2(6)),  &
       (SweepNum2,Hed2(7)), (Azm2,Hed2(8)), (Tilt2,Hed2(9)), (Alt2,Hed2(10)), (Rdel2,Hed2(11)), (Alt_Sel2,Hed2(12)), &
       (Ra2,Hed2(13))

  !  This is a extremely high powered interpolation subroutine that
  !  uses two rays at a time to objectively analyze radial data to a Cartesian
  !  grid.
  ! . . . . . . . . Written originally by Dave Jorgensen...12/1/83
  ! . . . . . . . . Updated for for/aft scanning 4/20/86
  ! . . . . . . . . Converted to HP 9000 (UNIX) and updated 6/24/93
  ! . . . . . . . . Updated to work with SOLO "mrd" files 5/5/97
  ! . . . . . . . . Added topography clutter removal 1/20/04
  
  !  Both reflectivity and Doppler velocity are interpolated to cartesian
  !  space.  X and Y are filled, Z is interpolated.

  !  Array Hed contains the header data
  !  Array Vel contains the radial velocity data
  !  Array Refl contains the reflectivity data
  
  !  First, read in the topo data (if available)

  If (Itopo .eq. 1) Then
     Nrec = 1
     
     Do j = 1, Max_topoy
        Read (66,Iostat=Ierr,Rec=Nrec) (Elev(i,j),i=1,Max_topox)
        Nrec = Nrec + 1
     End Do
     
     Close (66)
  End If
  
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
  
  Iyi = Init_Time(1) + Iadjyr
  Ini = Init_Time(2)
  Idi = Init_Time(3)
  Ihi = Init_Time(4)
  Imi = Init_Time(5)
  Isi = Init_Time(6)
  
  T_Init = Timz(Ihi,Imi,Isi)
  If (Idi .gt. Isday) T_Init = T_Init + 86400.0
  Stime  = Timz(Ihs,Ims,Iss)
  If (Ids .gt. Isday) Stime = Stime + 86400.0
  
1000 Call GetBeam (Itime1, Hed1, Vel1, Ref1, Ierr)

  If (Ierr .ne. 0) Then
     Write (LM,'(" Error returning from GetBeam:",i5)') Ierr
     Return
  End If
  
  !  Make sure the data is within the time window
  
  Iy = Itime1(1)
  In = Itime1(2)
  Id = Itime1(3)
  Ih = Itime1(4)
  Im = Itime1(5)
  Is = Itime1(6)
  
  If (     TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iys,Ins,Ids,Ihs,Ims,Iss)) Go To 1000

  If (.not.TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iye,Ine,Ide,Ihe,Ime,Ise)) Then
     Write (LM,'(" End time found:",1x,i2.2,"/",i2.2,"/",i4,1x,3i2.2/ &
          " Start Time:",i2.2,"/",i2.2,"/",i4,1x,3i2.2,/ &
          "   End Time:",i2.2,"/",i2.2,"/",i4,1x,3i2.2)') &
          In,   Id,  Iy,  Ih,  Im,  Is, Ins, Ids, Iys, Ihs, Ims, Iss, Ine, Ide, Iye, Ihe, Ime, Ise
    Return
 End If

 Time1 = Timz(Ih,Im,Is)
 If (Id .gt. Isday) Time1 = Time1 + 86400.0

 If (Azraw1 .lt. Sazm(1)) Go To 1000
 If (Azraw1 .gt. Eazm(1) .and. Azraw1 .gt. Eazm(2)) Go To 1000

 If (Nmosm .eq. 1) Then
    Olat1 = Olat + Sv*(Time1-T_Init)/ 111190.0
    Olon1 = Olon + Su*(Time1-T_Init)/(111190.0*Cos(Olat*Ctr))
    Tdif = Time1 - T_Init
 Else
    Olat1 = Olat
    Olon1 = Olon
    Tdif = Time1 - Stime
 End If

 Dlat = Plat1 - Olat1
 Dlon = Plon1 - Olon1
 Xz1= 111.19 * Cos(Plat1*Ctr) * Dlon
 Yz1= 111.19 * Dlat
 N_Gates1 = Zgat1
 Zz1 = Alt1/1000.0

 ! Avoid using gates contaminated by ground/sea return.
 ! Calculate the bin value corresponding to point where the bottom
 ! of beam drops below surface; assume half-power beamwidth of 1.0 deg.
 !  If aircraft height is expressed in terms of pressure altitude,
 ! use the digital terrain map (if available)

 If (Azraw1 .gt. 88.95 .and. Azraw1 .lt. 271.05 .and. Alt_Sel1 .eq. 0.0) Then
    Sn = Sign(1.0,Sin(Azm1*Ctr))
    Cs = Cos((Azm1+Sn) * Ctr)  *  Cos (Tilt1 * Ctr)
    A = -Alt1/Cs
    B = Rdel1 * 1000.0
    C = Gtln1 * 1000.0
    Ngatsfc = Nint((A - B)/C)

    If (Ngatsfc .le. N_Gates1) Then
       Do i = Ngatsfc, N_Gates1
          Vel1(i) = -999.0
          Ref1(i) = -999.0
       End Do
    End If
 End If

 !  Remove the sfc clutter with the aid of a digital topo map

 If (Azraw1 .gt. 88.95 .and. Azraw1 .lt. 271.05 .and. Alt_Sel1 .eq. 1.0 .and. Itopo .eq. 1) Then
    Call Terrain (Azm1, Ang1, Alt1, Ra1, N_Gates1, Rdel1, Gtln1, Plat1, Plon1, Elev, Vel1, Ref1)
 End If

 !  Print out some info at the start of each sweep

 Call Ctme (Time1,Ih,Im,Is)

 Write (LM,10) Ih, Im, Is, Plat1, Plon1, Xz1, Yz1, Alt1, Tilt1, Olat1, Olon1, SweepNum1, Rdel1, Tdif, Alt_Sel1, &
      Gtln1, Zgat1, Azm1

10 Format (/1X,3I2.2,' Lat:',F7.3,' Lon:',F8.3,' Xz:',F6.2,' Yz:', F6.2,' Alt:',F5.0,' Tilt:',F5.1/ &
        3X,' Olat:',F7.3,' Olon:',F8.3,' Sw #:',F4.0,' Rdel:', F6.2,' Del Tm:',F6.0/3x,' Alt S:',f2.0,' Gtln:',F6.3, &
        ' Gates:',F4.0,' Azm:',F6.1)

2000 Call GetBeam (Itime2, Hed2, Vel2, Ref2, Ierr)

 If (Ierr .ne. 0) Then
    Write (LM,'(" Error returning from GetBeam:",i5)') Ierr
    Return
 End If

 !  Make sure the data is within the time window

 Iy = Itime2(1)
 In = Itime2(2)
 Id = Itime2(3)
 Ih = Itime2(4)
 Im = Itime2(5)
 Is = Itime2(6)

 If (     TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iys,Ins,Ids,Ihs,Ims,Iss)) Go To 2000

 If (.not.TimeOnBefore(Iy,In,Id,Ih,Im,Is,Iye,Ine,Ide,Ihe,Ime,Ise)) Then

    Write (LM,'(" End time found:",1x,i2.2,"/",i2.2,"/",i4,1x, 3i2.2,1x," End Time:",i2.2,"/",i2.2,"/",i4,1x,3i2.2)') &
         In,Id,Iy,Ih,Im,Is, Ine,Ide,Iye,Ihe,Ime,Ise
    Return
 End If

 Time2 = Timz(Ih,Im,Is)
 If (Id .gt. Isday) Time2 = Time2 + 86400.0
 
 !  Data within the azimuth window?
 
 If (Azraw2 .gt. Eazm(1) .and. Sazm(2) .eq. 0.0) Go To 1000
 If (Azraw2 .gt. Eazm(1) .and. Azraw2 .lt. Sazm(2)) Go To 3000
 If (Azraw2 .gt. Eazm(1) .and. Azraw2 .gt. Eazm(2)) Go To 1000
 
 !  New sweep?
 
 If (SweepNum1 .ne. SweepNum2) Go To 1000
 
 !  Time jump?
 
 If (Time2-Time1 .gt. 2.0) Go To 3000

 !  Move the grid

 If (Nmosm .eq. 1) Then
    Olat2 = Olat + Sv*(Time2-T_Init)/ 111190.0
    Olon2 = Olon + Su*(Time2-T_Init)/(111190.0*Cos(Olat*Ctr))
    Tdif = Time2 - T_Init
 Else
    Olat2 = Olat
    Olon2 = Olon
    Tdif = Time2 - Stime
 End If
 
!  Compute the aircraft's position in the grid

 Dlat = Plat2 - Olat2
 Dlon = Plon2 - Olon2
 Xz2 = 111.19 * Cos(Plat2*Ctr) * Dlon
 Yz2 = 111.19 * Dlat
 Xz = (Xz2 + Xz1) * .5
 Yz = (Yz2 + Yz1) * .5
 
 CosAng1 = Cos(Ang1*Ctr)
 SinAng1 = Sin(Ang1*Ctr)
 CosAng2 = Cos(Ang2*Ctr)
 SinAng2 = Sin(Ang2*Ctr)
 
 SinAng1_Ang2 = SinAng1  - SinAng2
 CosAng1_Ang2 = CosAng1  - CosAng2
 
 CosAzm1 = Cos(Azm1*Ctr)
 SinAzm1 = Sin(Azm1*Ctr)
 CosAzm2 = Cos(Azm2*Ctr)
 SinAzm2 = Sin(Azm2*Ctr)
 
 Xval1 = SinAzm1 * SinAng1
 Xval2 = SinAzm2 * SinAng2
 Yval1 = SinAzm1 * CosAng1
 Yval2 = SinAzm2 * CosAng2
 Xval = (Xval1 + Xval2) / 2.0
 Yval = (Yval1 + Yval2) / 2.0
 Zval1 = CosAzm1
 Zval2 = CosAzm2
 
 N_Gates2 = Zgat2
 Zz2 = Alt2/1000.0
 C1 = (1-Zval1**2) * 5.885E-05
 C2 = (1-Zval2**2) * 5.885E-05
 
 ! Avoid using gates below ground/sea level.
 ! Calculate the bin value corresponding to point where the beam
 ! drops below surface; assume half beamwidth of 1.0 deg.
 
 !  If aircraft height is expressed in terms of pressure altitude,
 !  the digital topo data are used (if available)
 
 If (Azraw2 .gt. 88.95 .and. Azraw2 .lt. 271.05 .and. Alt_Sel1 .eq. 0.0) Then
    Sn = Sign(1.0, 180.0-Azm2)
    Cs = Cos((Azm2+Sn)*Ctr)
    A = -Alt2/Cs
    B = Rdel2 * 1000.0
    C = Gtln2 * 1000.0
    Ngatsfc = Nint((A - B)/C)
    
    If (Ngatsfc .le. N_Gates2) Then
       Do i = Max0(1,Ngatsfc), N_Gates2
          Vel2(i) = -999.0
          Ref2(i) = -999.0
       End Do
    End If
 End If
 
 !     Remove the rough terrain via a digital topo map

 If (Azraw2 .gt. 88.95 .and. Azraw2 .lt. 271.05 .and. Alt_Sel2 .eq. 1.0 .and. Itopo .eq. 1) Then
    Call Terrain (Azm2, Ang2, Alt2, Ra2, N_Gates2, Rdel2, Gtln2, Plat2, Plon2, Elev, Vel2, Ref2)
 End If
 
 N_Gates = Min0(N_gates1,N_Gates2)
 
 Do i = 1, N_Gates
    
    !  Check if the bin contains a valid velocity
    
    If (Vel1(i) .eq. -999.0) Cycle
    If (Vel2(i) .eq. -999.0) Cycle
    
    !  Compute the range to the bin
    
    R = Rdel2 + Float(i-1) * Gtln2

    !  Is the range within the specified limit?
    
    If (R .gt. RangeMax) Cycle
    
    !  Compute the y (North-South) location
    
    Y = Yz + Yval * R
    Jy = Nint(Y/Sy)
    If (Jy .lt. 1 .or. Jy .gt. Jmax) Cycle
    
    !  Compute the x (East-West) location
    
    X = Xz + Xval * R
    Ix = Nint(X/Sx)
    If (Ix .lt. 1 .or. Ix .gt. Imax) Cycle
    
    !  Compute the height of the bin
    
    Rsq = R*R
    Z1 = Zz1 + R*Zval1 + Rsq*C1
    Z2 = Zz2 + R*Zval2 + Rsq*C2
    Zbig = Amax1(Z1,Z2)
    Zsml = Amin1(Z1,Z2)
    
    !  Get the reflectivity (Z) and velocity (m/s) for the two bins
    
    Q1 = 10.0**(Ref1(i)/10.0)
    Q2 = 10.0**(Ref2(i)/10.0)
    V2 = Vel2(i)
    V1 = Vel1(i)
    
    If (Abs(Z1-Z2) .gt. .00001) Then
       Do Kz = 1, Kmax
          Hgt = Float(Kz-1)*Sz + Z0
          
          !  Require grid point to be bounded above and below by actual radial
          !  data, thus avoiding vertical extrapolation.
          
          If (Hgt .ge. Zsml .and. Hgt .le. Zbig) Then
             Diff  = (Hgt-Z2)/(Z1-Z2)
             Vradi = V2 + Diff*(V1-V2)
             Azmi  = Azm2 + Diff*(Azm1-Azm2)
             Tilti = Tilt2 + Diff*(Tilt1-Tilt2)
             Angi  = Atn2((SinAng2 + Diff*(SinAng1_Ang2)), (CosAng2 + Diff*(CosAng1_Ang2)))/Ctr
             Zi   = Q2 + Diff*(Q1-Q2)
             If (Zi .gt. -999.0)  Zi = 10.0*Alog10(Zi)

             Call PutData (Ix, Jy, Kz, Data, Angi, Azmi, Vradi, Zi, Tilti, R, Tdif, Itype)
          End If
       End Do
    End If
 End Do

3000 Do i = 1, N_Gates2
    Vel1(i) = Vel2(i)
    Ref1(i) = Ref2(i)
 End Do

 N_Gates1=N_Gates2
 Xz1 = Xz2
 Yz1 = Yz2
 Time1 = Time2
 
 Do i = 1, 10
    Hed1(i) = Hed2(i)
 End Do
 
Go To 2000
End Subroutine Cint

Subroutine GetBeam (Itime, Hed, Vel, Ref, Ierr)

  ! This Subroutine Reads the merged Doppler and Reflectivity data from
  ! disc (.cor file) that has been already unfolded and cleaned up
  
  !  Array HED contains the header data needed by CINT:
  
  !  Hed(1)  = Compass azimuth of beam relative to North
  !  Hed(2)  = Latitude [deg]
  !  Hed(3)  = Longitude [deg]
  !  Hed(4)  = Gate Length [km]
  !  Hed(5)  = # of gates
  !  Hed(6)  = original azimuth corrected for roll
  !  Hed(7)  = Sweep Number
  !  Hed(8)  = track relative rotation angle from zenith
  !  Hed(9)  = track relative tilt angle (fore/aft tilt)
  !  Hed(10) = altitude [m]
  !  Hed(11) = Range delay [km]             as AGL (using radar altitude)
  !  Hed(12) = 0 for RA, 1 for PA           or MSL (using pressure altitude)]
  !  Hed(13) = Original altitude when .cor was created (probably RA)
  
  !  One sweep consists of up to 400 radials of data with each radial
  !  containing up to 1024 gates.  The data is compacted by looking for
  !  the last good gate
  
  !  Format of the .cor file
  !   First a couple of header records:  
  !  The first record consists of a short header (45 integer*2 words) followed
  !  by a second header (40 integer*4 words)
  !  The data record then follows which consists of N_Gates of data.  The data
  !  thresholded prior to writing the file by marching inward from the end of the
  !  beam and finding the first good gate, calling that point N_Gates.  This
  !  methodology seems to work well for an airborne radar since the beam usually
  !  spends most of its time pointing into outer space or below the surface.
  
  !  The format of the disc header record is as follows:
  !   Ihead(1) = Hours
  !   Ihead(2) = Minutes
  !   Ihead(3) = Seconds
  !   Ihead(4) = Rotation Azimuth [deg*10] from straight up (not roll corrected)
  !   Ihead(5) = Latitude [deg]
  !   Ihead(6) = Latitude [min]
  !   Ihead(7) = Latitude [sec*10.0]
  !   Ihead(8) = Longitude [deg]
  !   Ihead(9) = Longitude [min]
  !   Ihead(10) = Longitude [sec*10.0]
  !   Ihead(11) = Altitude [m] defined at load time as AGL (i.e., radar
  !               altitude) or MSL (i.e., pressure altitude appropriate
  !               for work over highly variable terrain)
  !   Ihead(12) = Roll [deg*10.0]
  !   Ihead(13) = Head [deg*10.0]
  !   Ihead(14) = Drift [deg*10.0]
  !   Ihead(15) = Pitch [deg*10.0]
  !   Ihead(16) = Tilt angle [deg*10]  fore/aft pointing angle raw recorded value
  !   Ihead(17) = Nyquist velocity [m/s*10]
  !   Ihead(18) = Julian date
  !   Ihead(19) = # of azimuth samples
  !   Ihead(20) = Gate length [m]
  !   Ihead(21) = Range delay [m]
  !   Ihead(22) = Ground speed [m/s*64.0]
  !   Ihead(23) = Vertical airspeed [m/s*64.0]
  !   Ihead(24-27) = Flight number (8 ASCII characters)
  !   Ihead(28-33) = Storm name (12 ASCII characters)
  !   Ihead(34) = Wind dir [deg*10]
  !   Ihead(35) = Generation date
  !   Ihead(36) = Wind speed [m/s*10]
  !   Ihead(37) = Threshold [dBm] used to identify "noise"
  !   Ihead(38) = "corrected" tilt angle for pitch, roll, and drift
  !   Ihead(39) = # of good gates 
  !   Ihead(40) = radial velocity correction for the ray due to ground speed
  !   Ihead(41) = Sweep Number
  !   Ihead(42) = Maximum number of gates
  !   Ihead(43) = Tilt correction flag (0=no; 1=yes) if data has been changed
  !   Ihead(44) = flag for altitude 0 for RA, 1 for PA
  !   Ihead(45) = Not Used
  
  !   The second header currently uses only the first 6 words:
  !     Iheadr(1) = year (00 - 99)
  !     Iheadr(2) = month (00 - 12)
  !     Iheadr(3) = day (00 - 31)
  !     Iheadr(4) = hour (00 - 23)
  !     Iheadr(5) = minute (00 - 59)
  !     Iheadr(6) = second (00 - 59)

  !   The velocity and reflectivity data are scaled for each gate as
  !   a single 16 bit integer (to save disc space) as follows:

  !                                 bit #
  !                             1111110000000000
  !                             5432109876543210
  !                             zzzzzzvvvvvvvvvv
  !    where:
  !          z: scaled reflectivity [6 bits] 0 - 63 dBZ
  !          v: velocity data [10 bits] +-68 m/s 0.15 m/s resolution
  !       If signal <noise then word is =-1
  
  Common /Two/ Sx, Sy, Sz, Z0, Olat, Olon, Sazm(2), Eazm(2), Nmosm, Su, Sv, Xniq
  Common /Tre/ Flid, Project
  Common /Write/ LM

  Character*8 Flid
  Character*16 Project
  
  Dimension Hed(*), Vel(*), Ref(*), Itime(*)
  
  Integer*2 Ihead(45), Ibufr(1024)
  Integer*4 Iheadr(40)
  
  !  Read the header first
  
1 Read (75,End=1000,Err=999,Iostat=Ierr) (Ihead(i),i=1,45), (Iheadr(i),i=1,40)

  If (Ihead(41) .eq. -1) Go To 1    !   psuedo end of sweep mark
  
  N_Gates = Ihead(39)
  Max_Gates = Ihead(42)
  
  If (N_Gates .gt. 1024 .or. Max_Gates .gt. 1024 .or. N_Gates .gt. Max_Gates) Then
     Write (LM,'("Problem with n_gates or max_gates",2i9)') N_Gates, Max_Gates
     Stop
  End If

  Itime(1) = Iheadr(1)   ! year
  Itime(2) = Iheadr(2)   ! month
  Itime(3) = Iheadr(3)   ! day
  Itime(4) = Iheadr(4)   ! hour
  Itime(5) = Iheadr(5)   ! minute
  Itime(6) = Iheadr(6)   ! second
  
  
  !  Next read the data record
  !  Note that only the first N_Gates out of Max_Gates are actually stored,
  !  Since thresholding is performed prior to writing the file.
  
  Read (75,End=1000,Err=999,Iostat=Ierr) (Ibufr(i),i=1,N_Gates)
  
  !  Load the header information
  
  !  Latitude

  Hed(2) = Float(Ihead(5)) + Float(Ihead(6))/60.0 + Float(Ihead(7))/36000.0
  
  !  Longitude
  
  Hed(3) = Float(Ihead(8)) + Float(Ihead(9))/60.0 + Float(Ihead(10))/36000.0
  
  !  Gate Length
  
  Hed(4) = Float(Ihead(20))/1000.0
  
  !  Number of good gates in the ray
  
  Hed(5) = Float(Ihead(39))
  
  !  Sweep Number
  
  Hed(7) = Float(Ihead(41))
  
  !  Altitude (m)
  
  Hed(10) = Float(Ihead(11))
  
  !  Range delay (km) "RDEL"
  
  Hed(11) = Float(Ihead(21))/1000.0
  
  !  Which altitude flag
  
  Hed(12) = Float(Ihead(44))  ! =0 for RA, = 1 for PA
  
  !  Original altitude
  
  Hed(13) = Float(Ihead(45))
  
  !  Project and flight identifiers
  
  Write (Flid,'(4a2)') (Ihead(kk),kk=24,27)
  Write (Project,'(6a2,4H- - )') (Ihead(kk),kk=28,33)
  
  !  Nyquist interval (m/s)
  
  Xniq = Float(Ihead(17)) / 10.0
  
  !  Correct the ray for bad antenna behavior and compute angles
  !  relative to the track
  
  Call Correct_Ray (Vel, Ref, Ihead, Tilt, Rot, CompAz, Azm_RC)
  
  !  Azimuth roll corrected
  
  Hed(6) = Azm_RC
  
  !  Rotation angle from zenith
  
  Hed(8) = Rot
  
  !  Track relative tilt angle (fore/aft direction)
  
  Hed(9) = Tilt
  
  !  Compass azimuth from North
  
  Hed(1) = CompAz
  
  !  Load up the data
  
  Do i = 1, Max_Gates
     If (i .le. N_Gates) Then
        Ival = Ibufr(i)
     Else
        Ival = -1
     End If
     
     If (Ival .eq. -1) Then
        Vel(i) = -999.0
        Ref(i) = -999.0
     Else
        Ivel = Ibits(Ival,0,10)
        Vel(i) = Float(Ivel - 511) / 7.5
        Ref(i) = Ibits (Ival,10,6)
     End If
  End Do
  
  Return
  
999 Write (LM,'("Error:",i5)') Ierr
  Return
  
1000 Write (LM,'(/" End of file on input disc file",i5)') Ierr
  Return
  
End Subroutine GetBeam

Subroutine Correct_Ray (Vel, Ref, Ihead, Tilt_TR, Rot_TR, CompAz, Rot_RC)
  
  !  Subroutine to "correct" the ray for various effects
  
  !  Variables within the program are:
  !   Vel:   Array of radial velocity values
  !   Ref:   Array of reflectivity values
  !   Tilt:  Track relative tilt angle (fore/aft direction) + is forward
  !   Rot:   Rotation angle from zenith
  !   CompAz:Beam pointing angle relative to North (Compass Azimuth)
  
  Dimension Vel(*), Ref(*)
  
  Integer*2 Ihead(*)
  
  Common /Six/ Ctr, Eloff, CantDrft, Rotcr, PitchCor, GsCor, RangeMax 
  Common /Tre/ Flid, Project
  
  Character*8 Flid
  Character*16 Project
  
  !  Rotation angle [deg]
  
  Rot_Raw = Rotcr + Float(Ihead(4))/10.0
  
  !  Tilt angle [deg], raw (fore/aft direction)
  
  Tilt = Float(Ihead(16))/10.0 - 180.0
  
  !  Roll angle [deg]
  
  Roll = Float(Ihead(12))/10.0
  
  !  Correct the rotation angle for roll
  
  Rot_RC = Amod(Rot_Raw+Roll+360.0,360.0)
  
  !  Pitch angle [deg], corrected for radar mounting error
  
  Ptc = PitchCor + Float(Ihead(15))/10.0
  
  !  Gate Length [km]
  
  Gtln = Float(Ihead(20))/1000.0
  
  !  Range delay [km]
  
  Rdel = Float(Ihead(21))/1000.0
  
  !  Altitude [km]

  Alt = Float(Ihead(11))/1000.0

  !  Drift angle [deg]

  Dft = Float(Ihead(14))/10.0

  !  Heading [deg]

  Head = Float(Ihead(13))/10.0

  !  Track [deg]

  Trck = Amod(Head+Dft+360.0,360.0)

  !  Ground speed [m/s], corrected for INE bias/drift

  Gs = GsCor + Float(Ihead(22))/64.0

  !  Vertical aircraft speed [m/s]

  Vv = Float(Ihead(23))/64.0

  !  Maximum number of gates

  N_Gates = Ihead(39)

  Call Trans_Coords (Rot_RC, Tilt, Dft, Ptc, Trck, CompAz, Rot_TR, Tilt_TR)

  !  Compute the "correction" to the radial velocities caused by the
  !  aircraft's motion
  !  First make sure that this step was not done when the data were loaded
  !  on the disc
  
  Old_Cor = Float(IHead(40))/10.0
  Tilt_TR = Tilt_TR + Eloff - CantDrft*Sin(Rot_Raw*Ctr)
  
  If (Ihead(43) .eq. 0) Old_Cor = 0.0
  
  Cor = -Sin(Tilt_TR*Ctr)*Gs - Cos(Rot_TR*Ctr) * Vv
  Corvl = Cor - Old_Cor
  
  !     Recalculate the radial velocity correction in all but the ELDORA
  !     data
  
  If (Flid(7:7) .ne. 'E') Then
     Do i = 1, N_Gates
        If (Vel(i) .gt. -999.0) Then
           Vel(i) = Vel(i) + Corvl
        End If
     End Do
  End If
  
  Return
  
End Subroutine Correct_Ray

Subroutine PutData (i,j,k, Data, Ang, Azm, Vel, Z, Tilt, Range, Tdif, Itype)
  
  ! Ang (deg) = ground-based azimuth, i.e., compass angle for beam
  !     projected onto horizontal plane (0 to 359.999)
  ! Azm (deg) = azimuth in spherical coordinates
  ! Tilt (deg) = elev relative to track
  
  Common /One/ Imax, Jmax, Kmax
  Common /Six/ Ctr, Eloff, CantDrft, Rotcr, PitchCor, GsCor, RangeMax 

  Dimension Data(7,Imax,Jmax,Kmax)

  Character*3 Itype
  
  !  Subroutine to store data values
  !   Values falling within a grid volume are average if the angles are
  !   within 20 degrees
  
  !     Data(1,i,j,k) = compass angle  [degrees] relative to north
  !     Data(2,i,j,k) = radial velocity  [m/s]
  !     Data(3,i,j,k) = azimuth  [degrees] relative to zenith
  !     Data(4,i,j,k) = reflectivity [dBZ]
  !     Data(5,i,j,k) = counter
  !     Data(6,i,j,k) = Range (km) to bin
  !     Data(7,i,j,k) = Time relative to nominal
  
  !  Which direction is the antenna pointing? (fore or aft)
  
  !  store only the correct tilt (fore or aft)
  
  If (Tilt .lt. 0.0 .and. Itype .eq. 'for') Return
  If (Tilt .gt. 0.0 .and. Itype .eq. 'aft') Return
  
  
  !  Look at the data currently stored for the grid point
  
  Ang_Old  = Data(1,i,j,k)
  Vel_Old  = Data(2,i,j,k)
  Azm_Old  = Data(3,i,j,k)
  Ref_Old  = Data(4,i,j,k)
  Xnum     = Data(5,i,j,k)
  Rng_Old  = Data(6,i,j,k)
  Tdif_Old = Data(7,i,j,k)
  
  !  Insure that Azm and Ang are within 20 degrees of currently stored
  !  values to avoid averaging data from disparate directions
  
  If (Xnum .gt. 0.0) Then
     If (Abs(Azdif(Azm_Old,Azm,180.0)) .gt. 20.0) Return
     If (Abs(Azdif(Ang_Old,Ang,180.0)) .gt. 20.0) Return
  End If
  
  !  Compute the new range (km) to the bin
  
  R = (Rng_Old*Xnum + Range)/(Xnum+1.0)
  
  !  Compute the new velocity
  
  Velocity = (Vel_Old*Xnum + Vel)/(Xnum+1.0)
  
  !  Compute the new azimuth
  
  Azimuth = (Azm_Old*Xnum + Azm)/(Xnum+1.0)
  
  !  compute the new compass angle from north
  
  X = Cos(Ang_Old*Ctr)
  Y = Sin(Ang_Old*Ctr)
  X1= Cos(Ang*Ctr)
  Y1= Sin(Ang*Ctr)
  XT= (X*Xnum + X1)/(Xnum+1.0)
  YT= (Y*Xnum + Y1)/(Xnum+1.0)
  Angle = Atn2(Yt,Xt)/Ctr
  
  !  Compute the average time
  
  Tdif1 = (Tdif_Old*Xnum + Tdif)/(Xnum+1.0)
  
  !  Use the max reflectivity of all the gates in the bin
  
  Refz = Ref_Old
  If (Z .gt. Ref_Old) Refz = Z
  
  Xnum = Xnum + 1
  
  Data(1,i,j,k) = Angle
  Data(2,i,j,k) = Velocity
  Data(3,i,j,k) = Azimuth
  Data(4,i,j,k) = Refz
  Data(5,i,j,k) = Xnum
  Data(6,i,j,k) = R
  Data(7,i,j,k) = Tdif1
  
  Return
End Subroutine PutData

Subroutine GetData (i,j,k, Data, Ang, Vel, Azm, Z, Rng, Tdif)
  
  !  This Subroutine returns the velocity vector information
  !  The format of the array is as follows:

  !      angle from north = Data(1)
  !      velocity         = Data(2)
  !      azimuth angle    = Data(3)
  !      reflectivity     = Data(4)
  !      count            = Data(5)
  !      range            = Data(6)
  !      time difference  = Data(7)
  
  Common /One/ Imax, Jmax, Kmax
  
  Dimension Data(7,Imax,Jmax,Kmax)
  
  Ang = Data(1,i,j,k)
  Vel = Data(2,i,j,k)
  Azm = Data(3,i,j,k)
  Z   = Data(4,i,j,k)
  Xnum= Data(5,i,j,k) 
  Rng = Data(6,i,j,k)
  Tdif= Data(7,i,j,k)
  
  If (Xnum .eq. 0.0) Then
     Ang = -999.0
     Vel = -999.0
     Azm = -999.0
     Z   = -999.0
     Rng = -999.0
     Tdif= -999.0
     Return
  End If
  
  Return
End Subroutine GetData

Subroutine Trans_Coords (Rotation, Tilt_Raw, Drft, Ptch, Trck, Azimuth, Elev_Zen, Tilt_TR)
  
  !  Routine to translate the aircraft relative tail-radar measured angles
  !  to track relative angles
  
  !  Input parameters:
  !       Rotation:  Rotation angle of the antenna from zenith 
  !                  (used to be called "azimuth") corrected for roll
  !       Tilt_Raw:  Fore/Aft angle measured relative to a normal plane
  !                  to the airframe
  !       Drft:      Drift angle
  !       Ptch:      Pitch angle
  !       Trck:      Track angle
  
  !  Output Parameters:
  !       Azimuth:   Compass direction of the beam measured from North
  !       Elev_Zen:  Elevation angle of the beam measured from zenith
  !       Tilt_TR:   Track relative tilt angle (Fore/Aft)
  
  Common /Six/ Ctr, Eloff, CantDrft, Rotcr, PitchCor, GsCor, RangeMax 
  
  !  Convert angles to radians
  
  Rot = Rotation * Ctr
  Tlt = Tilt_Raw * Ctr
  Dri = Drft * Ctr
  Pit = Ptch * Ctr
  Trk = Trck * Ctr
  
  !  Direction cosine for x (distance normal to the track)
  
  x = cos(Rot)*sin(Dri)*sin(Pit)*cos(Tlt) + cos(Dri)*sin(Rot)*cos(Tlt) - sin(Dri)*cos(Pit)*sin(Tlt)
  
  !  Direction cosine for y (distance along track due to fore/aft pointing)
  
  y = -cos(Rot)*cos(Dri)*sin(Pit)*cos(Tlt) + sin(Dri)*sin(Rot)*cos(Tlt) + cos(Dri)*cos(Pit)*sin(Tlt)
  
  !  Direction cosine for z (height)
  
  z = cos(Pit)*cos(Rot)*cos(Tlt) + sin(Pit)*sin(Tlt)
  
  !  Track relative tilt angle (fore/aft looking angle)
  
  Tilt_TR = Asin(y)/Ctr
  
  !  Azimuth (x,y) angle (compass direction) relative to the track
  
  Azm_TR = Atn2(x,y)/Ctr
  
  !  Azimuth angle relative to North
  
  Azimuth = Amod(Azm_TR+Trck,360.0)
  
  !  Elevation angle from the horizontal (+ or - 90 deg)
  
  Elev_Hor = Asin(z)/Ctr
  
  !  Elevation angle from zenith (0 to 180 deg)
  
  Elev_Zen = 90.0 - Elev_Hor
  
  Return
End Subroutine Trans_Coords

Subroutine Terrain (Azm, Ang, Alt, Ra, Ngates, Rdel, Gtln, Olat, Olon, Topo, Vel, Ref)
  
  !     Routine to remove the gates below the ground
  !     To account for beamwidth, angles within 1 degree of the pointing
  !     angle are used
  
  ! Azm is the rotation angle from straight up
  ! Ang is the pointing angle relative to North
  ! Alt is the pressure altitude (m)
  ! Ra  is the radar altitude (m)
  ! Ngates is the number of the last good gate
  ! Rdel is the range delay (km)
  ! Gtln is the gate length (m)
  ! Olat, Olon is the postition of the aircraft
  ! Topo is the data array holding the topo info (m)
  ! Vel is the array containing the radial velocity data
  ! Ref is the array containng the reflectivity data
  ! Tlat, Tlon is the origin (lower left corner) of the topo array
  ! Tsx, Tsy is the grid spacing (km) of the top array
  !     Max_topox, Max_topoy is the number of grid points (x,y) in the
  !     topo array

  ! Thus the topo array covers Max_topox*Tsx, Max_topoy*Tsy (x,y) kmxkm
  
  Common /Topo/ Tlat, Tlon, Tsx, Tsy, Max_topox, Max_topoy, Itopo
  Common /Write/ LM
  
  Dimension Topo (Max_topox, Max_topoy)
  Dimension Vel(Ngates), Ref(Ngates)
  
  Azm1 = Azm - 1.0          ! One side of the beam
  Azm2 = Azm + 1.0          ! The other side of the beam (beamwidth is about 2 deg)
  
  !  First, find the topo height at the AC location
  
  AC_topo = Terrain_Hgt(180.0, 0.0, Alt, 256, Rdel, Gtln, Olat, Olon, Topo, Jgate, Gate_Hgt)
  
  !  Add the topo height at the AC location to the AC Ra to get height
  !  of the AC above MSL
  
  AC_Hgt = Ra + AC_topo
  
  !  Find the terrain height and gate number where the bottom of the
  !  beam hits the ground
  
  Thgt1 = Terrain_Hgt (Azm1, Ang, AC_Hgt, Ngates, Rdel, Gtln, Olat, Olon, Topo, Jgate1, Gate_Hgt1)
  
  !  Lastly, find the terrain height and gate number where the top of
  !  the beam hits the ground
  
  Thgt2 = Terrain_Hgt (Azm2, Ang, AC_Hgt, Ngates, Rdel, Gtln, Olat, Olon, Topo, Jgate2, Gate_Hgt2)
  
  !  Define the gate of the lowest part of the beam to intercept the ground
  
  NgateJ = JGate1
  If (Gate_Hgt2 .le. Gate_Hgt1) NgateJ = JGate2
  
  If (NgateJ .lt. NGates) Then
     Do i = NgateJ, NGates
        Vel(i) = -999.0
        Ref(i) = -999.0
     End Do
  End If
  
  Igates = Ngates - NgateJ + 1
  
  If (Igates .gt. 0) Then
     Write (LM,'("Removing ",i5," gates from Azm: ",F6.2," due to ", "terrain")') Igates, Azm
     !         Write (LM,*) Ang, AC_Hgt, Ngates, NgateJ, Gate_Hgt1, Gate_Hgt2, Jgate1, Jgate2, Thgt1, Thgt2
  End If
  
  Return
  
End Subroutine Terrain

Function Terrain_Hgt(Azm, Ang, Alt, Ngates, Rdel, Gtln, Olat, Olon, Topo, Jgate, Gate_Hgt)
  
  Parameter (Rearth = 6366.8056, Pi = 3.14159)
  Parameter (Ctr = Pi/180.0)
  
  Common /Topo/ Tlat, Tlon, Tsx, Tsy, Max_topox, Max_topoy, Itopo
  Common /Write/ LM
  
  Dimension Topo (Max_topox, Max_topoy)
  
  !     Function to return the height of the terrain at the point the beam
  !     intercepts the ground
  
  Sang = Sin(Ang*Ctr)
  Cang = Cos(Ang*Ctr)
  Sazm = Sin(Azm*Ctr)
  Cazm = Cos(Azm*Ctr)
  
  Xval = Sazm * Sang
  Yval = Sazm * Cang
  Zval = Cazm
  
  ZZ = Alt/1000.0
  C2 = (1.0 - (Zval*Zval)) * 5.885E-05
  Jgate = 9999
  Gate_Hgt = -9999

  Terrain_Hgt = -999.0  

  Do i = 1, NGates
     R = Rdel + Float(i-1)*Gtln
     X = Xval*R
     Y = Yval*R
     Z = ZZ + R*Zval + R*R*C2
     
     !     to get the lon/lat for the gate position (X,Y)
     
     Clat = Olat + (Asin(Y/Rearth)/Ctr)
     Clon = Olon + (Asin(X/(Rearth*Cos(Clat*Ctr)))/Ctr)
     
     !     Get the maxi terrain height from 9 grid points surrounding (X,Y)
     
     Dx = 111.19 * Cos(Tlat*Ctr) / Tsx * (Clon-Tlon)
     ip = Nint(Dx + 1.0)
     Dy = 111.19 / Tsy * (Clat-Tlat)
     jp = nint(Dy + 1.0)
     
     !  Is the gate location within the topo domain?
     
     If (Ip .lt. 2 .or. Ip .gt. Max_topox-1 .or. Jp .lt. 2 .or. Jp .gt. Max_topoy-1) Then
        Write (LM,'("Gate location is outside the topo domain", " (Ip,Jp): ",2i6," (Max_topox, Max_topoy): ",2I5)') &
             Ip, Jp, Max_topox, Max_topoy
        Return
     End If

     Theight = 0.0
     
     Do jj = jp-1, jp+1
        Do ii = ip-1, ip+1
           Terrain_Hgt = Amax1(Theight, Topo(ii,jj))
        End Do
     End Do
     
     !     Once the height of the gate is lower than the terrain height then exit
     
     Jgate = i
     Gate_Hgt = Z * 1000.0
     If (Z*1000.0 .le. Theight) Return
  End Do

  !     If this return is taken then the good data in the beam never hit
  !     the ground.  Bogus in a gate number greater than Ngates to
  !     distinguish this event
  
  Jgate = Ngates + 1
  
  Return
End Function Terrain_Hgt

Block Data
   
   Common /One/ Imax, Jmax, Kmax
   Common /Two/ Sx, Sy, Sz, Z0, Olat, Olon, Sazm(2), Eazm(2), Nmosm, Su, Sv, Xniq
   Common /Tre/ Flid, Project
   Common /Five/ Itime_Limits(12), Init_Time(6), Isday, Iadjyr
   Common /Six/ Ctr, Eloff, CantDrft, Rotcr, PitchCor, GsCor, RangeMax 
   Common /Topo/ Tlat, Tlon, Tsx, Tsy, Max_topox, Max_topoy, Itopo
   Common /Write/ LM
   
   Character*8 Flid
   Character*16 Project
   
   Data Flid /'????????'/, Project /'????????????????'/
   Data Ctr/.01745329/, Eloff /0.0/, CantDrft /0.0/, Rotcr /0.0/
   Data PitchCor /0.0/, GsCor /0.0/, LM /60/
   
End Block Data
