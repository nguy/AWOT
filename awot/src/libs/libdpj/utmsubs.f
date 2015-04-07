      subroutine constants

      common /consts/   PI, FOURTHPI, deg2rad, rad2deg
      common /ellipses/ ellipsnum(24), ellipsname(24), ell_eqrad(24),
     $     ell_eccsq(24)

      integer ellipsnum
      character*20 ellipsname
      real*8 ell_eqrad
      real*8 ell_eccsq
      
      real*8 PI
      real*8 FOURTHPI
      real*8 deg2rad
      real*8 rad2deg
      
      PI = 3.14159265
      FOURTHPI = PI / 4
      deg2rad = PI / 180
      rad2deg = 180.0 / PI
      
      data ellipsnum /-1,1,2,3,4,5,6,7,8,9,10,11,
     +     12,13,14,15,16,17,18,19,20,
     +     21,22,23/
      
      data ellipsname /"Placeholder",
     +     "Airy",
     +     "Australian National",
     +     "Bessel 1841",
     +     "Bessel 1841 (Nambia)",
     +     "Clarke 1866",
     +     "Clarke 1880",
     +     "Everest",
     +     "Fischer1960(Mercury)",
     +     "Fischer 1968",
     +     "GRS 1967",
     +     "GRS 1980",
     +     "Helmert 1906",
     +     "Hough",
     +     "International",
     +     "Krassovsky",
     +     "Modified Airy",
     +     "Modified Everest",
     +     "Modified Fischer1960",
     +     "South American 1969",
     +     "WGS 60",
     +     "WGS 66",
     +     "WGS-72",
     +     "WGS-84"/
      
      data ell_eqrad / 0.,
     +     6377563.,
     +     6378160.,
     +     6377397.,
     +     6377484.,
     +     6378206.,
     +     6378249.,
     +     6377276.,
     +     6378166.,
     +     6378150.,
     +     6378160.,
     +     6378137.,
     +     6378200.,
     +     6378270.,
     +     6378388.,
     +     6378245.,
     +     6377340.,
     +     6377304.,
     +     6378155.,
     +     6378160.,
     +     6378165.,
     +     6378145.,
     +     6378135.,
     +     6378137/

      data ell_eccsq / 0.,
     +     0.00667054,
     +     0.006694542,
     +     0.006674372,
     +     0.006674372,
     +     0.006768658,
     +     0.006803511,
     +     0.006637847,
     +     0.006693422,
     +     0.006693422,
     +     0.006694605,
     +     0.00669438,
     +     0.006693422,
     +     0.00672267,
     +     0.00672267,
     +     0.006693422,
     +     0.00667054,
     +     0.006637847,
     +     0.006693422,
     +     0.006694542,
     +     0.006693422,
     +     0.006694542,
     +     0.006694318,
     +     0.00669438/

      return
      end

      subroutine lltoutm (ReferenceEllipsoid, Lat, Long, 
     +     UTMNorthing, UTMEasting, UTMZone, ZoneNumber)

      common /consts/   PI, FOURTHPI, deg2rad, rad2deg
      common /ellipses/ ellipsnum(24), ellipsname(24), ell_eqrad(24),
     $     ell_eccsq(24)

      integer ellipsnum
      character*20 ellipsname
      real*8 ell_eqrad
      real*8 ell_eccsq
      
      real*8 PI
      real*8 FOURTHPI
      real*8 deg2rad
      real*8 rad2deg
      
      integer ReferenceEllipsoid
      real*8 Lat
      real*8 Long
      real*8 UTMNorthing
      real*8 UTMEasting
      character UTMZone
      character*1 UTMLetterDesignator

c-----converts lat/long to UTM coords.  Equations from USGS Bulletin 1532 
c-----East Longitudes are positive, West longitudes are negative. 
c-----North latitudes are positive, South latitudes are negative
c-----Lat and Long are in decimal degrees

      real*8 a
      real*8 eccSquared
      real*8 k0
      
      real*8 LongOrigin
      real*8 eccPrimeSquared
      real*8 N, T, C, BA, M
      real*8 LatRad
      real*8 LongRad
      real*8 LongTemp
      real*8 LongOriginRad
      integer ZoneNumber

c---Initialize the constants

      Call Constants

      a = ell_eqrad(ReferenceEllipsoid)
      eccSquared = ell_eccsq(ReferenceEllipsoid)
      k0 = 0.9996

c-----Make sure the longitude is between -180.00 .. 179.9

      longt = (Long+180.0)/360.0
      LongTemp = (Long+180.0)-longt*360.0-180.0

c----- -180.00 .. 179.9

      LatRad = Lat*deg2rad
      LongRad = LongTemp*deg2rad

      longt = (LongTemp + 180.0)/6.0
      ZoneNumber = longt + 1
  
      if(Lat.ge.56.0.and.Lat.lt.64.0.and.LongTemp.ge.3.0.and.LongTemp.lt
     $     .12.0) ZoneNumber = 32

c----- Special zones for Svalbard

      if(Lat.ge.72.0 .and. Lat.lt.84.0 ) then
         if(LongTemp.ge.0.0  .and. LongTemp.lt. 9.0 ) ZoneNumber = 31
         if(LongTemp.ge.9.0  .and. LongTemp.lt.21.0 ) ZoneNumber = 33
         if(LongTemp.ge.21.0 .and. LongTemp.lt.33.0 ) ZoneNumber = 35
         if(LongTemp.ge.33.0 .and. LongTemp.lt.42.0 ) ZoneNumber = 37
      end if
      
c-----+3 puts origin in middle of zone

      LongOrigin = (ZoneNumber - 1)*6 - 180 + 3
      LongOriginRad = LongOrigin*deg2rad

c-----compute the UTM Zone from the latitude and longitude

      UTMZone = UTMLetterDesignator(Lat)

      eccPrimeSquared = (eccSquared)/(1-eccSquared)
      
      N = a/sqrt(1-eccSquared*dsin(LatRad)*dsin(LatRad))
      T = tan(LatRad)*tan(LatRad)
      C = eccPrimeSquared*dcos(LatRad)*dcos(LatRad)
      BA = dcos(LatRad)*(LongRad-LongOriginRad)
      
      M = a*((1 - eccSquared/4 - 3*eccSquared*eccSquared/64 
     +     - 5*eccSquared*eccSquared*eccSquared/256)*LatRad
     +     - (3*eccSquared/8 + 3*eccSquared*eccSquared/32 
     +     + 45*eccSquared*eccSquared*eccSquared/1024)*
     +     dsin(2.d0*LatRad)
     +     + (15*eccSquared*eccSquared/256 
     +     + 45*eccSquared*eccSquared*eccSquared/1024)*
     +     dsin(4.d0*LatRad) 
     +     - (35*eccSquared*eccSquared*eccSquared/3072)*
     +     dsin(6.d0*LatRad))
      
      UTMEasting = k0*N*(BA+(1-T+C)*BA*BA*BA/6
     +     + (5-18*T+T*T+72*C-58*eccPrimeSquared)
     +     *BA*BA*BA*BA*BA/120)
     +     + 500000.0
      
      UTMNorthing = k0*(M+N*tan(LatRad)*(BA*BA/2+
     +     (5-T+9*C+4*C*C)*BA*BA*BA*BA/24
     +     + (61-58*T+T*T+600*C-330*eccPrimeSquared)
     +     *BA*BA*BA*BA*BA*BA/720))
      
      if (Lat .lt. 0.0) then

c-----10000000 meter offset for southern hemisphere

         UTMNorthing = UTMNorthing + 10000000.0
      end if
      
      return
      end

      function UTMLetterDesignator(Lat)
      real*8 Lat
      character*1 UTMLetterDesignator

c-----This routine determines the correct UTM letter designator for the given latitude
c-----returns 'Z' if latitude is outside the UTM limits of 84N to 80S

      if((84.ge.Lat) .and. (Lat.ge.72)) then
         UTMLetterDesignator = 'X'
      else if((72.gt.Lat) .and. (Lat.ge.64)) then
         UTMLetterDesignator = 'W'
      else if((64.gt.Lat) .and. (Lat.ge.56)) then
         UTMLetterDesignator = 'V'
      else if((56.gt.Lat) .and. (Lat.ge.48)) then
         UTMLetterDesignator = 'U'
      else if((48.gt.Lat) .and. (Lat.ge.40)) then
         UTMLetterDesignator = 'T'
      else if((40.gt.Lat) .and. (Lat.ge.32)) then
         UTMLetterDesignator = 'S'
      else if((32.gt.Lat) .and. (Lat.ge.24)) then
         UTMLetterDesignator = 'R'
      else if((24.gt.Lat) .and. (Lat.ge.16)) then
         UTMLetterDesignator = 'Q'
      else if((16.gt.Lat) .and. (Lat.ge.8)) then
         UTMLetterDesignator = 'P'
      else if(( 8.gt.Lat) .and. (Lat.ge.0)) then
         UTMLetterDesignator = 'N'
      else if(( 0.gt.Lat) .and. (Lat.ge.-8)) then
         UTMLetterDesignator = 'M'
      else if((-8.gt.Lat) .and. (Lat.ge.-16)) then
         UTMLetterDesignator = 'L'
      else if((-16.gt.Lat) .and. (Lat.ge.-24)) then
         UTMLetterDesignator = 'K'
      else if((-24.gt.Lat) .and. (Lat.ge.-32)) then
         UTMLetterDesignator = 'J'
      else if((-32.gt.Lat) .and. (Lat.ge.-40)) then
         UTMLetterDesignator = 'H'
      else if((-40.gt.Lat) .and. (Lat.ge.-48)) then
         UTMLetterDesignator = 'G'
      else if((-48.gt.Lat) .and. (Lat.ge.-56)) then
         UTMLetterDesignator = 'F'
      else if((-56.gt.Lat) .and. (Lat.ge.-64)) then
         UTMLetterDesignator = 'E'
      else if((-64.gt.Lat) .and. (Lat.ge.-72)) then
         UTMLetterDesignator = 'D'
      else if((-72.gt.Lat) .and. (Lat.ge.-80)) then
         UTMLetterDesignator = 'C'
      else 

c-----This is here as an error flag to show that the Latitude is outside the UTM limits
         UTMLetterDesignator = 'Z'
      end if


      return 
      end

      subroutine utmtoll (ReferenceEllipsoid, UTMNorthing, UTMEasting,
     $     UTMZone,Lat,  Long, ZoneNumber)

      common /consts/   PI, FOURTHPI, deg2rad, rad2deg
      common /ellipses/ ellipsnum(24), ellipsname(24), ell_eqrad(24),
     $     ell_eccsq(24)

      integer ellipsnum
      character*20 ellipsname
      real*8 ell_eqrad
      real*8 ell_eccsq
      
      real*8 PI
      real*8 FOURTHPI
      real*8 deg2rad
      real*8 rad2deg
      
      integer ReferenceEllipsoid
      real*8 UTMNorthing
      real*8 UTMEasting
      character UTMZone
      real*8 Lat
      real*8 Long
      integer ZoneNumber
      
      real*8 N1, T1, C1, R1, D, M
      real*8 LongOrigin
      real*8 mu, phi1, phi1Rad
      real*8 x, y
      character ZoneLetter

c-----1 for northern hemispher, 0 for southern

      integer NorthernHemisphere 
      real*8 k0
      real*8 a
      real*8 eccSquared
      real*8 eccPrimeSquared
      real*8 e1

c-----converts UTM coords to lat/long.  Equations from USGS Bulletin 1532 
c-----East Longitudes are positive, West longitudes are negative. 
c-----North latitudes are positive, South latitudes are negative
c-----Lat and Long are in decimal degrees. 

      Call Constants

      k0 = 0.9996
      a = ell_eqrad(ReferenceEllipsoid)
      eccSquared = ell_eccsq(ReferenceEllipsoid)
      e1 = (1-sqrt(1-eccSquared))/(1+sqrt(1-eccSquared))

c-----remove 500,000 meter offset for longitude

      x = UTMEasting - 500000.0 
      y = UTMNorthing

      ZoneLetter = UTMZone

      if((ichar(ZoneLetter) - ichar('N')).ge.0) then

c-----point is in northern hemisphere

         NorthernHemisphere = 1
      else

c-----point is in southern hemisphere
         NorthernHemisphere = 0

c-----remove 10,000,000 meter offset used for southern hemisphere

         y = y - 10000000.0
      end if

c-----+3 puts origin in middle of zone
      LongOrigin = (ZoneNumber - 1)*6 - 180 + 3
      
      eccPrimeSquared = (eccSquared)/(1-eccSquared)
      
      M = y / k0
      mu = M/(a*(1-eccSquared/4-3*eccSquared*
     +     eccSquared/64-5*eccSquared*
     +     eccSquared*eccSquared/256))
      
      phi1Rad = mu + (3*e1/2-27*e1*e1*e1/32)*
     +     dsin(2*mu) + (21*e1*e1/16-55*
     +     e1*e1*e1*e1/32)*dsin(4*mu)
     +     +(151*e1*e1*e1/96)*dsin(6*mu)
      
      phi1 = phi1Rad*rad2deg
      
      N1 = a/sqrt(1-eccSquared*dsin(phi1Rad)*
     +     dsin(phi1Rad))
      T1 = tan(phi1Rad)*tan(phi1Rad)
      C1 = eccPrimeSquared*dcos(phi1Rad)*dcos(phi1Rad)
      R1 = a*(1-eccSquared)/((1-eccSquared*
     +     dsin(phi1Rad)*dsin(phi1Rad))**1.5)
      D = x/(N1*k0)
      
      Lat = phi1Rad - (N1*tan(phi1Rad)/R1)*
     +     (D*D/2-(5+3*T1+10*C1-4*C1*C1-9*
     +     eccPrimeSquared)*D*D*D*D/24
     +     +(61+90*T1+298*C1+45*T1*T1-252*
     +     eccPrimeSquared-3*C1*C1)*D*D*D*D*D*D/720)
      Lat = Lat*rad2deg
      
      Long = (D-(1+2*T1+C1)*D*D*D/6+
     +     (5-2*C1+28*T1-3*C1*C1+
     +     8*eccPrimeSquared+24*T1*T1)
     +     *D*D*D*D*D/120)/dcos(phi1Rad)
      
      Long = LongOrigin + Long*rad2deg
      
      return
      end
