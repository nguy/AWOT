      !    "<23-Jun-1994 22:20:46><UTC>"
      Integer*4 Function Ichk(x,y,Xmin, Xmax, Ymin, Ymax)
! check to see if point is in box with boundaries xmin,xmax,ymin,ymax
! return 0: if in box
!        add 1 if less than xmin,
!        add 2 if greater than xmax,
!        add 4 if less than ymin,
!        add 8 if greater than ymax,
      Implicit none
      Real*4 x,y, Xmin, Xmax, Ymin, Ymax
  
      Ichk = 0
  
      If (x .lt. Xmin) Ichk = 1
      If (x .gt. Xmax) Ichk = 2
      If (y .gt. Ymax) Ichk = Ichk + 8
      If (y .lt. Ymin) Ichk = Ichk + 4
  
      Return
      End ! function Ichk ends
************************************************************************
