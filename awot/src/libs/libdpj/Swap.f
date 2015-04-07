      !    "<23-Jun-1994 22:20:48><UTC>"
      Subroutine Swap(x,y)
      Implicit none
! parameters 
      Real*4 x,y

! local variables
      Real*4 t
  
      t = x
      x = y
      y = t
  
      Return
      End ! subroutine Swap ends
************************************************************************
