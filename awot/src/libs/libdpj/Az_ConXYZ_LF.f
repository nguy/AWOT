C  "<08-Oct-1993 21:54:25><UTC>"
******************************************************************
      Subroutine Az_ConXYZ_LF (Azm, El, Track, Drift, Pitch, Roll,
     #      GrAz,GrEl,Xval, Yval, Zval)
  
c  Routine to translate the aircraft relative Lower Fuselage radar (LF)
c  measured angles to ground-based relative angles & direction cosines.
c  All angles passed from/to calling module are in degrees.
  
c  Input parameters:
c       Azm:  Rotation angle of the LF-antenna beam (called "azimuth").
c             Zero degrees is toward the nose, 90 is to the right wing.
c             Azm varies from 0 to 359.9999
c       El:   Angle measured relative to a plane containing the "horizontal"
c             axes of the airframe.  Plus is toward the ceiling.
c       Drift: Drift angle, angle from heading to track.  E.g., if
c             heading is 355, and track is 5, then drift is 10 degrees.
c       Pitch: Pitch angle, angle horizontal plane to P3's axis.  Normally
c             pitch is 2-3 degrees.
c       Roll: Roll angle, angle of rotation of P3 about  P3's axis.  Normally
c             roll is within a couple degrees of zero, but in sharp turns 
c             its magnitude can be near 45 degrees.  Left wing up is positive.
c       Track: Track angle, compass angle of P3's track.  Track varies from
c             0 (North) to 359.9999.
  
c  Output Parameters:
c       GrAz:  Compass direction of the beam measured from North, 
c             "ground-based" azimuth, 0 to 359.999 degrees.
c       GrEl:  Elevation angle of the beam measured,
c             "ground-based" elevation, +90 (zenith) to -90 (nadir) degrees.
c             This is the same as LF tilt!!
c       Xval,Yval,Zval: direction cosines of beam relative to earth
c             (track included).
  
      Implicit None
      Real*4 CTR  ! To convert between degrees and radians
      Parameter  (CTR = 3.1415926535897932/180.0)
      Real*4 Azm,El,Drift,Track,Pitch,Roll,GrAz,GrEl,Xval,Yval,Zval
      Real*4 AzmR,ElR,DriftR,PitchR,RollR,Angle,X,Y,Z
  
c  Convert angles to radians
  
      AzmR = Azm      * Ctr
      ElR = El       * Ctr
      DriftR = Drift * Ctr
      PitchR = Pitch * Ctr
      RollR = Roll * Ctr
  
  
c  Calculation is based on converting the direction cosines relative to
c  the P3's coordiate system to direction cosines relative to earth, i.e., 
c  ground-based reference.  This initial discussion will assume heading is
c  zero.  When the conversions have been finished, adjustments for heading=
c  track-drift will be made.  These comments "assume" the trignometric 
c  functions use degrees instead of radians.
c  I.  Let X0, Y0, Z0 be the direction cosines of the LF beam relative to the
c      P3's coordinate system.
c      Then X0 = Cos(El) * Sin(Azm)
c           Y0 = Cos(El) * Cos(Azm)
c           Z0 = Sin(El)
c      So X0**2 + Y0**2 + Z0**2 = 1
c  II. Direction cosines of the X,Y,Z axes of the P3 compensated for roll
c          with respect to the axes of the P3 not compensated for roll.
c          "Compensated for roll" is used to describe the P3 rotated about
c          its major axis so that wings are equal distances to the earth.
c      Let LX,MX,NX be the direction cosines of the X-axis 
c          P3 compensated for roll with respect to the P3's axes not 
c          not compensated for roll.
c      Let LY,MY,NY be the direction cosines of the Y-axis
c          P3 compensated for roll with respect to the P3's axes not 
c          not compensated for roll.
c      Let LZ,MZ,NZ be the direction cosines of the Z-axis
c          P3 compensated for roll with respect to the P3's axes not 
c          not compensated for roll.
c      Then LX = Cos(roll)
c           MX = 0
c           NX = Sin(roll)
c      Then LY = 0
c           MY = 1
c           NY = 0
c      Then LZ = -Sin(roll)
c           MZ = 0
c           NZ = Cos(roll) 
c  III. Calculation of beam's direction cosines in P3 compensated
c          for roll coordinates:
c      X1 = LX*X0 + MX*Y0 + NX*Z0
c         = Cos(roll)*Cos(El)*Sin(Azm)
c          +0+Sin(Roll)*Sin(El)
c      Y1 = LY*X0 + MY*Y0 + NY*Z0
c         =0 + Cos(El)*Cos(Azm) + 0
c      Z1 = LZ*X0 + MZ*Y0 + NZ*Z0
c         = - Sin(Roll)*Cos(El)*Sin(Az)
c          +0+Cos(Roll)*Sin(El)
c  IV. Direction cosines of the X,Y,Z axes (earth) with respect to the P3
c          axes (adjusted for roll!).
c      Let LambdaX,MuX,NuX be the direction cosines of the X-axis (earth) 
c          with respect to the P3's axes.
c      Let LambdaY,MuY,NuY be the direction cosines of the Y-axis (earth) 
c          with respect to the P3's axes.
c      Let LambdaZ,MuZ,NuZ be the direction cosines of the Z-axis (earth) 
c          with respect to the P3's axes.
c      Then LamdaX = 1
c           MuX = 0
c           NuX = 0
c      Then LamdaY = 0
c           MuY = Cos(Pitch) 
c           NuY = -Sin(Pitch)
c      Then LamdaZ = 0
c           MuZ = Sin(Pitch)
c           NuZ = Cos(Pitch)
c  V. Calculation of beam's direction cosines in earth coordinates:
c      X = LamdaX*X1 + MuX*Y1 + NuX*Z1
c        =Cos(roll)*Cos(El)*Sin(Azm)+Sin(Roll)*Sin(El) 
c        =Sin(Azm)*Cos(El)*Cos(Roll)
c         +Sin(El)*Sin(Roll)
c      Y = LamdaY*X1 + MuY*Y1 + NuY*Z1
c        = 0
c         +Cos(Pitch)**Cos(El)*Cos(Azm)
c         -Sin(Pitch)*(Cos(Roll)*Sin(El) -Sin(Roll)*Cos(El)*Sin(Azm))
c        =Cos(Azm)*Cos(El)*Cos(Pitch)
c         -Sin(El)*Sin(Pitch)*Cos(Roll)
c         +Sin(Azm)*Cos(El)*Sin(Pitch)*Sin(Roll)
c      Z = LamdaZ*X1 + MuZ*Y1 + NuZ*Z1
c        = 0 + Sin(Pitch)*Cos(El)*Cos(Azm) 
c         + Cos(Pitch)*(Cos(Roll)*Sin(El) -Sin(Roll)*Cos(El)*Sin(Azm))
c        =Cos(Azm)*Cos(El)*Sin(Pitch)
c         + Sin(El)*Cos(Pitch)*Cos(Roll)
c         - Sin(Azm)*Cos(El)*Cos(Pitch)*Sin(Roll)
  
  
      x = Sin(AzmR)*Cos(ElR)*Cos(RollR) +Sin(ElR)*Sin(RollR)
  
      y = Cos(AzmR)*Cos(ElR)*Cos(PitchR)
     >   -Sin(ElR)*Sin(PitchR)*Cos(RollR)
     >   +Sin(AzmR)*Cos(ElR)*Sin(PitchR)*Sin(RollR)
  
      z = Cos(AzmR)*Cos(ElR)*Sin(PitchR)
     >  + Sin(ElR)*Cos(PitchR)*Cos(RollR)
     >  - Sin(AzmR)*Cos(ElR)*Cos(PitchR)*Sin(RollR)
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
  
  
c  Azimuth (x,y) angle (compass direction) relative to the heading
      
      Angle = Atan2(x,y)/Ctr
  
c  Azimuth angle relative to North, heading included: "ground azimuth"
c  as heading = Track-Drift  
      GrAz     = Amod(Angle+Track-Drift+360.0,360.0)
  
c  Elevation angle from the horizontal (+ or - 90 deg): "ground-based elevation"
c  This is the same as LF Tilt.
  
      GrEl = Asin(z)/Ctr
      If (GrEl .gt. 90.0) GrEl = 90.0
      If (GrEl .lt. -90.0) GrEl = -90.0
  
  
      Xval = Cos(GrEl*Ctr) * Sin(GrAz*Ctr)
      Yval = Cos(GrEl*Ctr) * Cos(GrAz*Ctr)
      Zval = Z
      Return
      End ! Az_ConXYZ_LF ends
