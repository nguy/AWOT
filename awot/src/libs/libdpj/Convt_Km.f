      !    "<23-Jun-1994 22:20:52><UTC>"
      Subroutine Convt_Km(Plat,Plon,xkm,ykm,Clat,Clon)
  
c  Routine to convert Plat,Plon to kilometers on a flat Earth,
c      with Clat,Clon used as origin (i.e., xkm,ykm = 0.0,0.0).
  
      Implicit none
! parameters 
      Real*4 Plat,Plon,xkm,ykm,Clat,Clon

! local variables
      Real*4 Dlat,Dlon
  
      Dlon = Plon - Clon
      Dlat = Plat - Clat
      xkm = Dlon * 111.19 * Cos(Plat * 0.01745329)
      ykm = Dlat * 111.19
  
      Return
      End  ! subroutine Convt_km ends
************************************************************************
