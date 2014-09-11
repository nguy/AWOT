C  "<08-Oct-1993 21:54:34><UTC>"
************************************************************************
      Real*4 Function Tilt_LF (Azm, El, Pitch, Roll)
  
c  Routine to calculate Lower Fuselage radar (LF) tilt.
c  All angles passed from/to calling module are in degrees.
  
c  Input parameters:
c       Azm:  Rotation angle of the LF-antenna beam (called "azimuth").
c             Zero degrees is toward the nose, 90 is to the right wing.
c             Azm varies from 0 to 359.9999
c       El:   Angle measured relative to a plane containing the "horizontal"
c             axes of the airframe.  Plus is toward the ceiling.
c       Pitch: Pitch angle, angle horizontal plane to P3's axis.  Normally
c             pitch is 2-3 degrees.
c       Roll: Roll angle, angle of rotation of P3 about  P3's axis.  Normally
c             roll is within a couple degrees of zero, but in sharp turns 
c             its magnitude can be near 45 degrees.  Left wing up is positive.
  
c  Output:
c       Tilt_LF:  Elevation angle of the beam measured,
c             "ground-based" elevation, +90 (zenith) to -90 (nadir) degrees.
  
      Implicit None
      Real*4 CTR  ! To convert between degrees and radians
      Parameter  (CTR = 3.1415926535897932/180.0)
      Real*4 Azm,El,Pitch,Roll,Z
      Real*4 AzmR,ElR,PitchR,RollR
  
c  Convert angles to radians
  
      AzmR = Azm      * Ctr
      ElR = El       * Ctr
      PitchR = Pitch * Ctr
      RollR = Roll * Ctr
  
  
c  Please refer to the comments in Az_ConXYZ_LF for derivation of
c  formula to calculate z, the direction cosine of the beam in earth
c  coordinates (z is up).
  
      z = Cos(AzmR)*Cos(ElR)*Sin(PitchR)
     >  + Sin(ElR)*Cos(PitchR)*Cos(RollR)
     >  - Sin(AzmR)*Cos(ElR)*Cos(PitchR)*Sin(RollR)
  
      if (z .gt. 1.0) then
	  z = 1.0
      else if (z .lt. -1.0) then
	  z = -1.0
      end if
  
  
c  Elevation angle from the horizontal (+ or - 90 deg): "ground-based elevation"
c  This is the same as LF Tilt.
  
      Tilt_LF = Asin(z)/Ctr

      Return
      end ! Function Tilt_LF ends
