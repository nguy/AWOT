      !    "<23-Jun-1994 22:20:45><UTC>"
 
      Subroutine Clipit(Xkm1,Xkm2,Ykm1,Ykm2,Xmin, Xmax, Ymin, Ymax
     >   ,ScX,ScY)
      
! Plot line segment after it has been clipped! 
      Implicit none

      Real*4 Xkm1,Xkm2,Ykm1,Ykm2,
     >    Xmin, Xmax, Ymin, Ymax, ScX, ScY

! local variables
      Real*4  x,x1,x2,y,y1,y2

      Integer*4 Ic1,Ic2,Ichk
 
      x1 = Xkm1
      x2 = Xkm2
      y1 = Ykm1
      y2 = Ykm2
 
   10 Ic1 = Ichk(x1,y1,Xmin, Xmax, Ymin, Ymax)
      Ic2 = Ichk(x2,y2,Xmin, Xmax, Ymin, Ymax)
 
C  Check if both points are outside window
 
      If (Iand(Ic1,Ic2) .ne. 0) Return
  
C  Check if both points are inside window
  
      If (Ic1 .eq. 0 .and. Ic2 .eq. 0)Then
          x = x2 * ScX
          y = y2 * ScY
          Call Plot (x,y,3)
          x = x1 * ScX
          y = y1 * ScY 
          Call Plot (x,y,2)
          Return
      End If
  
C  If point #1 is inside, swap points
  
      If (Ic1 .eq. 0) Then
        Call Iswap(Ic1,Ic2)
        Call Swap(x1,x2)
        Call Swap(y1,y2)
      End If
  
C  Did boundary cross x (longitude) border?
  
      If (Iand(Ic1,3) .ne. 0) Go To 1
  
C  Went through latitude, move outside point to latitude
  
      y = Ymax
      If (Iand(Ic1,12) .eq. 4) y = Ymin
      x1 = x1 + (y-y1)*(x2-x1)/(y2-y1)
      y1 = y
  
      Go To 10
  
C  Went through longitude, move outside point to longitude
  
    1 continue
  
      x = Xmin
      If (Iand(Ic1,3) .eq. 2) x = Xmax
      y1 = y1 + (x-x1)*(y2-y1)/(x2-x1)
      x1 = x
 
      Go To 10
  
      End  ! subroutine clipit ends
***********************************************************************
