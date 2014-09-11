      !    "<23-Jun-1994 22:20:47><UTC>"
      Subroutine Iswap(i,j)
      Implicit none
! parameters 
      Integer*4 i,j

! local variables
      Integer*4 it
  
      it = i
      i = j
      j = it
  
      Return
      End ! subroutine Iswap ends
************************************************************************
