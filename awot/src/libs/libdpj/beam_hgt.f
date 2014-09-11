      Function Beam_Hgt(Range, Elev)

c     Routine to calculate the radar beam height using the 4/3 Earth
C     radius technique

c     Input variables:
c     Range - Range from radar [km]
c     Elev - Elevation angle [deg] from horizontal

c     Beam_Hgt is returned as [km]

      Data Rah/5.885E-5/, Ctr/0.0174532925/

c  The elevation angle term

      Term1 = Range * Sin(Elev*Ctr)

c  The 4/3 Earth radius of curvature term

      Term2 = Range * Range * Rah * Cos(Elev*Ctr) * Cos(Elev*Ctr)
      Beam_Hgt = Term1 + Term2

      Return
      End
