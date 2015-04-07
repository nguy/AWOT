      Subroutine Pltsg

C Plot One Segment - Straight Line, Dashed Line, Or Label

      Common /Clotx/ Es, Fs, Ek, Fk, Spi, Spj, Kount, Ipen, Klab, Cvar,
     . Dnum, Idsh, Pn, Qn, Imumd, Mumd, Numd, Hgt
  
c  This subroutine was modified by R. Hueftle Jan 1988 so that there
c     would be a more uniform clearance around the label, i.e., the
c     line segments are close but not intersect: this will at least hold
c     true for points plotted until a point (Fs,Es) is
c     outside a circle containing the label.
c   (Xc,Yc) is the center of the circle and Radc is the radius of the label
c     circle (all in terms of Plot coordinates).
  
      Parameter (Pi=3.1415926535897)                                   !RAH
      Parameter (CharW=.75) ! approx. ratio of character cell width to Hgt!RAH
      Parameter (Clear=.80) ! ratio of minimum label clearance to Hgt  !RAH
 
      Save Xc, Yc, Radc                                                !RAH

      Hgt2=Hgt/2.0
      Fsspj = Fs*Spj
      Esspi = Es*Spi
      Fsfk = Fs - Fk
      Ang = Atan2(Spi*(Es - Ek), Spj*Fsfk)
      If (Kount .Ne. 10) Go To 1
      If (Ipen .Eq. 3) Go To 1
      If (Klab .Eq. 0) Go To 1
      If (Idsh .Eq. 1 .And. Klab .Eq. 3) Go To 1
      If (Idsh .Eq. 0 .And. Klab .Eq. 2) Go To 1
      Ipen = 3  
                                                                        !RAH
c Get number of characters in label to determine size of label circle.  !RAH
                                                                        !RAH
      Nchar=1                                                           !RAH
      If (Cvar.lt.0) Nchar = 2                                          !RAH
      If (Abs(Cvar).Ge. 9.5)                                            !RAH
     1      Nchar = Nchar + Int ( Alog10(Abs(Cvar)+.5) )                !RAH
                                                                        !RAH
      Radc = Hgt*(Clear+Float(Nchar)*CharW*.5)  ! radius of clearance   !RAH
      Fkspj = Fk*Spj                                                    !RAH
      Ekspi = Ek*Spi                                                    !RAH
  
c  The location for the start of the label needs to be made.  Assume    !RAH
c     that the character cell width = CharW * Hgt.  Start label         !RAH
c     Hgt*Clear from last point with pen down,                          !RAH
c     i.e., (Fkspj,Ekspi).  (Xn,Yn) is the bottom-left corner of the    !RAH
c     label.  The function Atan2 returns a value in range -Pi to Pi.    !RAH
c     To keep the labels heads-up some adjustments need to be made for  !RAH
c     some of the values from  Atan2.                                   !RAH
  
      If ((Ang.Lt.-Pi/2. ).or.(Ang.Gt.Pi/2.)) then                      !RAH
          Ang = Ang + Pi ! move from 3->9 o'clock to 9->3 o'clock       !RAH
          C = Cos (Ang)                                                 !RAH
          S = Sin (Ang)                                                 !RAH
                                                                        !RAH
c  Position Xn,Yn (for room for full label) for bottom corner.          !RAH
                                                                        !RAH
          Xn = Fkspj - C*Hgt*(Clear+Float(Nchar)*CharW) + S*Hgt2        !RAH
          Yn = Ekspi - S*Hgt*(Clear+Float(Nchar)*CharW) - C*Hgt2        !RAH
                                                                        !RAH
c     Xc, Yc is center of circle of radius Radc containing the label    !RAH
                                                                        !RAH
          Xc = Fkspj - C*Radc                                           !RAH
          Yc = Ekspi - S*Radc                                           !RAH
      Else                                                              !RAH
          C = Cos(Ang)                                                  !RAH
          S = Sin(Ang)                                                  !RAH
                                                                        !RAH
c     Xn,Yn will be the bottom corner point to start the label (Cvar).  !RAH
                                                                        !RAH
          Xn=Fkspj + C*Hgt  + S*Hgt2                                    !RAH
          Yn=Ekspi + S*Hgt  - C*Hgt2                                    !RAH
                                                                        !RAH
c     Xc, Yc is center of circle of radius Radc containing the label    !RAH
                                                                        !RAH
          Xc = Fkspj + C*Radc                                           !RAH
          Yc = Ekspi + S*Radc                                           !RAH
      End If                                                            !RAH
      An = Ang * 180./Pi                                                !RAH
      Call Numbr(Xn,Yn,Hgt   ,Cvar,An ,-1)  ! plot label                !RAH
 
 1    If (Idsh .Eq. 0) Call Plot(Fsspj, Esspi, Ipen)
      If (Idsh .Eq. 1) Call Dash (Fsspj, Esspi, Ipen, Pn, Qn)
      If (Imumd .Eq. 0) Go To 3
      Mumd = Numd
      If (Fsfk .Ne. 0.) Mumd = Float(Numd)/(1. +
     . .6366197723*Ang*(Spi/Spj - 1.))
 3    Fk = Fs
      Ek = Es
      Kount = Kount + 1
***   If (Kount .Gt. Mumd + 10) Ipen = 2                                !RAH
  
c   Turn on pen if there has been label and current point is outside.   !RAH
  
      If (Kount .Gt. 10) Then                                           !RAH
          If (Sqrt( (Fsspj-Xc)**2 + (Esspi-Yc)**2) .Ge. Radc)           !RAH
     1        Ipen = 2                                                  !RAH
      End If                                                            !RAH
  
      Return
      End
