C  "<08-Oct-1993 21:54:26><UTC>"
****************************************************************************
      Subroutine Az_Tilt_ConXYZ (Azm, El, Track, Drift, Pitch,
     #      GrAz,GrEl,Tilt,Xval, Yval, Zval)
  
c  Routine to translate the aircraft relative tail-radar measured angles
c  to ground-based relative angles & direction cosines.
c  All angles passed from/to calling module are in degrees.
  
c  Input parameters:
c       Azm:  Rotation angle of the tail-antenna beam (called "azimuth"),
c             corrected for roll.  E.g., a raw azimuth of 350 and a roll of 
c             15 would give an Azm of 5 degrees. Azm varies from 0 to 359.9999.
c       El:  Fore/Aft angle measured relative to a normal plane
c             to the airframe, (called "raw elevation").  
c             This value is positive for fore, negative for aft.  In
c             P3 tapes before July88 the raw elevation was recorded 180 plus
c             our current raw elevation .  In P3 flights
c             before 1989 El was usually within a few degrees of zero
c             in attempt to keep Tilt very close to zero.  Now El often is
c             +/- 25 for the FAST dual doppler scanning techniques.
c       Drift: Drift angle, angle from heading to track.  E.g., if
c             heading is 355, and track is 5, then drift is 10 degrees.
c       Pitch: Pitch angle, angle horizontal plane to P3's axis.  Normally
c             pitch is 2-3 degrees.
c       Track: Track angle, compass angle of P3's track.  Track varies from
c             0 (North) to 359.9999.
  
c  Output Parameters:
c       GrAz:  Compass direction of the beam measured from North, 
c             "ground-based" azimuth, 0 to 359.999 degrees.
c       GrEl:  Elevation angle of the beam measured,
c             "ground-based" elevation, +90 (zenith) to -90 (nadir) degrees.
c       Tilt:   Track relative tilt angle (Fore +/Aft -).  Before 1989
c             tilt was kept close to zero as possible (to reduce
c             contamination of P3's ground speed in tail-radar 
c             radial velocity measurements), but with FAST scanning 
c             techniques tilt is usually +/- 25 degrees.
c       Xval,Yval,Zval: direction cosines of beam relative to earth
c             (track included).
  
       Implicit None
       Real*4 CTR  ! To convert between degrees and radians
       Parameter  (CTR = 3.1415926535897932/180.)
       Real*4 Azm,El,Drift,Track,Pitch,Tilt, GrAz,GrEl,Xval,Yval,Zval
       Real*4 AzmR,ElR,DriftR,PitchR,Angle,X,Y,Z
  
c  Convert angles to radians
  
      AzmR = Azm      * Ctr
      ElR = El       * Ctr
      DriftR = Drift * Ctr
      PitchR = Pitch * Ctr
  
  
c  Calculation is based on converting the direction cosines relative to
c  the P3's coordiate system to direction cosines relative to earth, i.e., 
c  ground-based reference.  This initial discussion will assume track is
c  zero.  When the conversions have been finished, adjustments for track
c  will be made.  These comments "assume" the trignometric functions
c  use degrees instead of radians.
c  I.  Let X0, Y0, Z0 be the direction cosines of the beam relative to the
c      P3's coordinate system.
c      Then X0 = Cos(El) * Sin(Azm)
c           Y0 = Sin(El)
c           Z0 = Cos(El) * Cos(Azm)
c      So X0**2 + Y0**2 + Z0**2 = 1
c  II. Direction cosines of the X,Y,Z axes (earth) with respect to the P3.
c      Let LambdaX,MuX,NuX be the direction cosines of the X-axis (earth) 
c          with respect to the P3's axes.
c      Let LambdaY,MuY,NuY be the direction cosines of the Y-axis (earth) 
c          with respect to the P3's axes.
c      Let LambdaZ,MuZ,NuZ be the direction cosines of the Z-axis (earth) 
c          with respect to the P3's axes.
c      Then LamdaX = Cos(drift)
c           MuX = Cos(Pitch) * Sin(-drift)
c           NuX = Sin(Pitch) * Sin(drift)
c      Then LamdaY = Sin(drift)
c           MuY = Cos(Pitch) * Cos(drift)
c           NuY = -Sin(Pitch) * Cos(drift)
c      Then LamdaZ = 0
c           MuZ = Sin(Pitch)
c           NuZ = Cos(Pitch)
c  III. Calculation of beam's direction cosines in earth coordinates:
c      X = LamdaX*X0 + MuX*Y0 + NuX*Z0
c        = Cos(Drift)*Cos(El)*Sin(Azm) 
c         -Cos(Pitch)*Sin(Drift)*Sin(El)
c         +Sin(Pitch)*Sin(Drift)*Cos(El)*Cos(Azm)
c      Y = LamdaY*X0 + MuY*Y0 + NuY*Z0
c        = Sin(Drift)*Cos(El)*Sin(Az)
c         +Cos(Pitch)*Cos(Drift)*Sin(El)
c         -Sin(Pitch)*Cos(Drift)*Cos(El)*Cos(Azm)
c      Z = LamdaZ*X0 + MuZ*Y0 + NuZ*Z0
c        = 0 + Sin(Pitch)*Sin(El) + Cos(Pitch)*Cos(El)*Cos(Azm)
c        
  
  
      x = Cos(DriftR)*Cos(ElR)*Sin(AzmR)
     >   -Cos(PitchR)*Sin(DriftR)*Sin(ElR)
     >   +Sin(PitchR)*Sin(DriftR)*Cos(ElR)*Cos(AzmR)
  
      y = Sin(DriftR)*Cos(ElR)*Sin(AzmR)
     >   +Cos(PitchR)*Cos(DriftR)*Sin(ElR)
     >   -Sin(PitchR)*Cos(DriftR)*Cos(ElR)*Cos(AzmR)
  
      z = Sin(PitchR)*Sin(ElR) + Cos(PitchR)*Cos(ElR)*Cos(AzmR)
  
      if (x .gt. 1.0) then
	  x = 1.0
      else if (x .lt. -1.0) then
	  x = -1.0
      end if
      if (y .gt. 1.0) then
	  y = 1.0
      else if (y .lt. -1.0) then
	  y = -1.0
      end if
      if (z .gt. 1.0) then
	  z = 1.0
      else if (z .lt. -1.0) then
	  z = -1.0
      end if

c  Tilt: angle to vertical plane that is perpendicular to track
  
      Tilt      = Asin(y)/Ctr
  
c  Azimuth (x,y) angle (compass direction) relative to the track
      
      Angle = Atan2(x,y)/Ctr
  
c  Azimuth angle relative to North, track included: "ground azimuth"
  
      GrAz     = Amod(Angle+Track+360.0,360.0)
  
c  Elevation angle from the horizontal (+ or - 90 deg): "ground-based elevation"
  
      GrEl = Asin(z)/Ctr
      If (GrEl .gt. 90.0) GrEl = 90.0
      If (GrEl .lt. -90.0) GrEl = -90.0
  
  
      Xval = Cos(GrEl*Ctr) * Sin(GrAz*Ctr)
      Yval = Cos(GrEl*Ctr) * Cos(GrAz*Ctr)
      Zval = Z
      Return
      End ! Az_Tilt_conXYZ ends
